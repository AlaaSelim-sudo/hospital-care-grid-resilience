"""Solve and replay a synchronized three-hospital surge during delayed restoration."""

import argparse
import json
from pathlib import Path

from run_fault_study import Bounds, LinearConstraint, Model, P, evaluate, linear_data, milp, np, read, scenario, write


POLICIES = ('autonomous', 'voltage_only', 'energy_first', 'coordinated')
EVENTS = ('dec2023', 'apr2024')
DIRECTORY = P / 'compound_case_data'


def design(event, pattern):
    restoration = dict(repair=21, tie_delay=12) if pattern in ('extended', 'reserve') else dict(repair=18, tie_delay=9)
    settings = scenario(event, **restoration)
    if pattern == 'reserve':
        settings['fuel_limit_l'] = (190., 190., 160.)
    fault = settings['fault_hour']
    original = settings['A'].copy()
    window = original[:, :, fault:fault + 12]
    cumulative = (np.cumsum(window, axis=2) + 1) // 2
    extra = np.diff(cumulative, axis=2, prepend=0)
    offsets = (0, 8, 16) if pattern == 'staggered' else (0, 0, 0)
    for h, offset in enumerate(offsets):
        start = fault + offset
        assert start + extra.shape[2] <= original.shape[2]
        settings['A'][h, :, start:start + extra.shape[2]] += extra[h]
    return settings, original, offsets


def objectives(model):
    essential = np.zeros(model.n)
    essential[model.names['plant']] = -1
    urgent, energy, timely = (np.zeros(model.n) for _ in range(3))
    for index, _, kind, arrival, start, duration in model.jobs:
        urgent[index] = -(kind == 0 and start - arrival <= 2)
        timely[index] = -(start - arrival <= (2, 8)[kind])
        energy[index] = -(duration * (4 if kind == 0 else 1.5) + (10 if kind == 0 else 0))
    # Final tie-break uses the same operating-cost terms as the main model.
    cost = np.zeros(model.n)
    for t in range(len(model.s['stages'])):
        hour = (model.s['start_hour'] + t) % 24
        cost[model.names['p'][:, t]] = .4 if 16 <= hour < 19 else (.12 if hour < 7 or hour >= 22 else .2)
    for name, price in {'generator': .36, 'generator_on': 1.8, 'charge': .02, 'discharge': .02}.items():
        cost[model.names[name]] = price
    return dict(essential=essential, urgent=urgent, energy=energy, timely=timely, cost=cost)


def solve_scipy(model, seconds, log_path):
    criteria = objectives(model)
    order = ('energy', 'urgent') if model.policy == 'energy_first' else ('urgent', 'energy')
    bounds = Bounds(model.lb, model.ub)
    constraints = [LinearConstraint(model.M, model.lo, model.hi)]
    logs, solution = [], None
    for name in ('essential', *order, 'timely', 'cost'):
        objective = criteria[name]
        result = milp(objective, integrality=model.integer, bounds=bounds, constraints=constraints,
                      options={'time_limit': seconds, 'mip_rel_gap': 0.002 if name == 'cost' else 0})
        log = model.log(result, name)
        if result.x is not None:
            solution = result.x
            if name != 'cost':
                value = round(float(objective @ solution))
                bound = log['dual_bound']
                upper = int(np.floor(-bound + 1e-5)) if bound is not None else None
                log.update(value=-value, upper_bound=upper, certified=upper == -value)
                constraints.append(LinearConstraint(objective, value - 1e-5, value + 1e-5))
        logs.append(log)
        log_path.write_text(json.dumps(logs, indent=2), encoding='utf-8')
        print(name, log, flush=True)
        if result.x is None:
            break
    return solution, logs


def solve(model, seconds, log_path, seed, stages=None, prime=False, criteria=None):
    import highspy as hs

    if criteria is None:
        criteria = objectives(model)
    order = ('energy', 'urgent') if model.policy == 'energy_first' else ('urgent', 'energy')
    solver = hs.Highs()
    solver.setOptionValue('output_flag', False)
    solver.setOptionValue('threads', 1)
    solver.setOptionValue('small_matrix_value', 1e-12)
    lp = hs.HighsLp()
    lp.num_col_, lp.num_row_ = model.n, model.M.shape[0]
    lp.col_cost_ = criteria['essential']
    lp.col_lower_, lp.col_upper_ = model.lb, model.ub
    lp.row_lower_, lp.row_upper_ = model.lo, model.hi
    lp.integrality_ = [hs.HighsVarType.kInteger if is_integer else hs.HighsVarType.kContinuous
                       for is_integer in model.integer]
    lp.a_matrix_.format_ = hs.MatrixFormat.kColwise
    lp.a_matrix_.start_, lp.a_matrix_.index_, lp.a_matrix_.value_ = model.M.indptr, model.M.indices, model.M.data
    assert solver.passModel(lp) == hs.HighsStatus.kOk
    columns = np.arange(model.n, dtype=np.int32)
    logs = []
    solution = np.load(seed)['x'] if seed else None
    if solution is not None:
        assert solution.shape == (model.n,)
        constraint_values = model.M @ solution
        assert np.all(constraint_values >= np.asarray(model.lo) - 1e-5)
        assert np.all(constraint_values <= np.asarray(model.hi) + 1e-5)
    for name in stages or ('essential', *order, 'timely', 'cost'):
        objective = criteria[name]
        solver.changeColsCost(model.n, columns, objective)
        solver.setOptionValue('mip_rel_gap', 0.002 if name == 'cost' else 0.)
        if prime and solution is not None:
            modes = np.concatenate([model.names[key].ravel() for key in ('plant', 'generator_on', 'battery_mode')]).astype(np.int32)
            fixed = np.rint(solution[modes])
            solver.changeColsBounds(len(modes), modes, fixed, fixed)
            solver.setOptionValue('time_limit', 12.)
            solver.run()
            hint = solver.getSolution()
            if hint.value_valid and objective @ hint.col_value < objective @ solution:
                solution = np.asarray(hint.col_value)
            solver.changeColsBounds(len(modes), modes, np.asarray(model.lb)[modes], np.asarray(model.ub)[modes])
        if solution is not None:
            assert solver.setSolution(model.n, columns, solution) == hs.HighsStatus.kOk
        budget = min(seconds, {'timely': 25, 'cost': 15}.get(name, seconds))
        elapsed = solver.getRunTime()
        solver.setOptionValue('time_limit', float(budget))
        solver.run()
        info = solver.getInfo()
        answer = solver.getSolution()
        status = solver.getModelStatus()
        log = dict(stage=name, status=int(status), message=solver.modelStatusToString(status),
                   solver=solver.version(), time_limit_s=budget, runtime_s=solver.getRunTime() - elapsed,
                   objective=float(info.objective_function_value) if answer.value_valid else None,
                   gap=float(info.mip_gap) if np.isfinite(info.mip_gap) else None,
                   dual_bound=float(info.mip_dual_bound) if np.isfinite(info.mip_dual_bound) else None)
        if seed:
            log['initial_seed'] = seed.name
        if answer.value_valid:
            candidate = np.asarray(answer.col_value)
            if solution is None or objective @ candidate <= objective @ solution + 1e-6:
                solution = candidate
            if name not in ('cost', 'diesel'):
                value = round(float(objective @ solution))
                bound = log['dual_bound']
                upper = int(np.floor(-bound + 1e-5)) if bound is not None else None
                log.update(value=-value, upper_bound=upper, certified=upper == -value)
                indices = np.flatnonzero(objective).astype(np.int32)
                assert solver.addRow(value, value, len(indices), indices, objective[indices]) == hs.HighsStatus.kOk
        logs.append(log)
        log_path.write_text(json.dumps(logs, indent=2), encoding='utf-8')
        print(name, log, flush=True)
        if not answer.value_valid:
            break
    return solution, logs


def record(event, pattern, policy, seconds, replay, seed=None, solver='highspy'):
    settings, original, offsets = design(event, pattern)
    guarded = policy in ('coordinated', 'energy_first')
    model = Model(settings, policy, linear_data(settings), tighten_v=.0003 if guarded else 0,
                  tighten_i=.012 if guarded else 0)
    stem = DIRECTORY / f'{event}_{pattern}_{policy}'
    log_path = stem.with_suffix('.json')
    if replay:
        logs = json.loads(log_path.read_text(encoding='utf-8'))
        solution_path = stem.with_suffix('.npz')
        solution = np.load(solution_path)['x'] if solution_path.exists() else None
    elif solver == 'scipy':
        solution, logs = solve_scipy(model, seconds, log_path)
    else:
        solution, logs = solve(model, seconds, log_path, seed)
    row = dict(event=event, pattern=pattern, policy=policy, baseline_arrivals=int(original.sum()),
               added_arrivals=int(settings['A'].sum() - original.sum()), fault_hour=settings['fault_hour'],
               isolation_h=settings['tie_delay_h'], repair_h=settings['repair_h'], xfm_kva=settings['xfm_kva'],
               tie_rating_A=settings['tie_rating_A'], surge_offsets_h='/'.join(map(str, offsets)))
    for log in logs:
        name = log['stage']
        for key in ('status', 'value', 'upper_bound', 'certified', 'gap'):
            row[f'{name}_{key}'] = log.get(key)
    if solution is not None:
        metrics, trace, _, currents, starts, _ = evaluate(settings, model, solution)
        assert metrics['linear_residual'] < 1e-4
        integers = solution[np.asarray(model.integer, dtype=bool)]
        assert np.max(np.abs(integers - np.rint(integers))) < 1e-4
        assert min(solution - model.lb) >= -1e-4 and min(model.ub - solution) >= -1e-4
        for name, vector in objectives(model).items():
            if name == 'cost':
                continue
            row[f'final_{name}'] = round(float(-vector @ solution))
            if row.get(f'{name}_value') is not None:
                assert row[f'final_{name}'] == row[f'{name}_value']
        row.update(metrics)
        prefix = stem.relative_to(P).as_posix()
        if replay:
            for suffix, reproduced in [('_replay.csv', trace), ('_current.csv', currents), ('_starts.csv', starts)]:
                archived_rows = read(prefix + suffix)
                assert len(archived_rows) == len(reproduced), (stem, suffix)
                for name, first in reproduced[0].items():
                    saved_column = [entry[name] for entry in archived_rows]
                    actual_column = [entry[name] for entry in reproduced]
                    if isinstance(first, (float, np.floating)):
                        assert np.allclose(np.asarray(saved_column, dtype=float), actual_column,
                                           atol=1e-6, rtol=1e-8), (stem, suffix, name)
                    else:
                        assert saved_column == list(map(str, actual_column)), (stem, suffix, name)
        else:
            np.savez_compressed(stem.with_suffix('.npz'), x=solution)
            write(prefix + '_replay.csv', trace)
            write(prefix + '_starts.csv', starts)
            write(prefix + '_current.csv', currents)
    if replay:
        archived = stem.with_suffix('.result.json')
        saved = json.loads(archived.read_text(encoding='utf-8'))
        assert saved.keys() == row.keys()
        for key, value in saved.items():
            if isinstance(value, float):
                assert np.isclose(value, row[key], atol=1e-6, rtol=1e-8), (stem, key, value, row[key])
            else:
                assert value == row[key], (stem, key, value, row[key])
    else:
        stem.with_suffix('.result.json').write_text(json.dumps(row, indent=2), encoding='utf-8')
    print(json.dumps(row), flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--event', choices=EVENTS)
    parser.add_argument('--policy', choices=POLICIES)
    parser.add_argument('--pattern', choices=('synchronized', 'staggered', 'extended', 'reserve'), default='synchronized')
    parser.add_argument('--seconds', type=int, default=70, help='Time limit per MILP objective')
    parser.add_argument('--seed', type=Path, help='Feasible schedule for the identical model constraints')
    parser.add_argument('--solver', choices=('highspy', 'scipy'), default='highspy')
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    if args.seed and args.solver == 'scipy':
        parser.error('--seed requires --solver highspy')
    DIRECTORY.mkdir(exist_ok=True)
    if args.replay:
        for event in EVENTS:
            simultaneous, original, _ = design(event, 'synchronized')
            staggered, _, _ = design(event, 'staggered')
            extended, _, _ = design(event, 'extended')
            reserve, _, _ = design(event, 'reserve')
            assert np.array_equal(simultaneous['A'].sum(axis=2), staggered['A'].sum(axis=2))
            assert np.array_equal(simultaneous['A'], extended['A'])
            assert np.array_equal(extended['A'], reserve['A'])
            assert np.all(extended['A'] >= original)
        rows = []
        for path in sorted(DIRECTORY.glob('*.result.json')):
            item = json.loads(path.read_text(encoding='utf-8'))
            rows.append(record(item['event'], item['pattern'], item['policy'], args.seconds, True))
        assert rows, 'No archived compound-stress results were found'
        write('Compound_Stress_Results.csv', rows)
    else:
        for event in (args.event,) if args.event else EVENTS:
            for policy in (args.policy,) if args.policy else POLICIES:
                record(event, args.pattern, policy, args.seconds, False, args.seed, args.solver)


if __name__ == '__main__':
    main()
