"""Long-isolation stress cases and a physical bound on essential-service loss."""

import argparse
import json
import math

from scipy.sparse import csr_matrix, vstack

from compound_stress import EVENTS, objectives, solve
from essential_reserve import size_reserve
from run_fault_study import Model, Network, P, background, evaluate, linear_data, np, read, scenario, write


DIRECTORY = P / 'extended_case_data'
# The mid-duration pair isolates the surge effect. The remaining cases vary
# duration and then initial stock, with the doubled surge held constant.
CASES = {'mid50': (15, .5, 160), 'mid100': (15, 1., 160),
         'long100': (18, 1., 160), 'long_reserve': (18, 1., 320),
         'long_no_surge': (18, 0., 320),
         'long_full': (18, 1., 160), 'long_full_reserve': (18, 1., 320)}
POLICIES = {
    'mid50': ('coordinated',), 'mid100': ('coordinated',),
    'long100': ('coordinated', 'voltage_only', 'energy_first'),
    'long_reserve': ('coordinated', 'voltage_only', 'energy_first'),
    'long_no_surge': ('voltage_only',),
    'long_full': ('coordinated',),
    'long_full_reserve': ('coordinated',),
}


def settings_for(event, name):
    isolation, surge, stock = CASES[name]
    settings = scenario(event, repair=isolation + 9, tie_delay=isolation)
    settings['fuel_limit_l'] = (stock, stock, 160)
    settings['backfeed_background_fraction'] = 1. if name.startswith('long_full') else .95
    fault = settings['fault_hour']
    window = settings['A'][:, :, fault:fault + 12]
    extra = np.diff(np.floor(surge * np.cumsum(window, axis=2) + .5).astype(int), axis=2, prepend=0)
    window += extra
    return settings, int(extra.sum())


def island_energy_bound(isolation, stock):
    # Maximize generator output over the integer number of hours committed.
    # The battery contributes at most its full usable discharge; demand,
    # losses beyond the discharge efficiency, and charging consume energy.
    generator = max(min(100 * hours, (stock - 1.5 * hours) / .30)
                    for hours in range(isolation + 1) if stock >= 1.5 * hours)
    return generator + .95 * (200 - 20)


def add_island_cut(model, isolation, stock):
    slots = min(isolation, math.floor(island_energy_bound(isolation, stock) / 60 + 1e-9))
    fault = model.s['fault_hour']
    cut = np.zeros((2, model.n))
    for h in range(2):
        cut[h, model.names['plant'][h, fault:fault + isolation]] = 1
    model.M = vstack([model.M, csr_matrix(cut)], format='csc')
    model.lo.extend([-np.inf, -np.inf])
    model.hi.extend([slots, slots])
    return 2 * (isolation - slots) * 60


def run(event, name, policy, seconds, replay, solved=None):
    settings, added = settings_for(event, name)
    isolation, surge, stock = CASES[name]
    guarded = policy != 'voltage_only'
    model = Model(settings, policy, linear_data(settings), tighten_v=.0003 if guarded else 0.,
                  tighten_i=.012 if guarded else 0.)
    physical_floor = add_island_cut(model, isolation, stock)
    stem = DIRECTORY / f'{event}_{name}_{policy}'
    log_path = stem.with_suffix('.json')
    vector_path = stem.with_suffix('.npz')
    if replay:
        logs = json.loads(log_path.read_text(encoding='utf-8'))
        solution = np.load(vector_path)['x'] if vector_path.exists() else None
    elif solved is not None:
        solution, logs = solved
    else:
        stages = ('essential', 'energy', 'urgent') if policy == 'energy_first' else ('essential', 'urgent')
        solution, logs = solve(model, seconds, log_path, seed=None, stages=stages, prime=True)
    pickup = settings.get('backfeed_background_fraction', 1.)
    deferred = (1 - pickup) * sum(p for _, p, _ in Network().base) * .25 * sum(
        background(settings['start_hour'] + t + k * .25)
        for t, stage in enumerate(settings['stages']) if stage == 'backfeed' for k in range(4))
    row = dict(event=event, case=name, policy=policy, isolation_h=isolation, repair_h=settings['repair_h'],
               surge_pct=int(surge * 100), extra_arrivals=added, fuel_per_island_l=stock,
               physical_deficit_floor_kwh=physical_floor, backfeed_native_fraction=pickup,
               deferred_native_kwh=deferred)
    for log in logs:
        for key in ('value', 'upper_bound', 'certified', 'gap', 'status'):
            row[f"{log['stage']}_{key}"] = log.get(key)
    if solution is not None:
        metrics, trace, _, currents, starts, _ = evaluate(settings, model, solution)
        assert metrics['linear_residual'] < 1e-4
        integer_values = solution[np.asarray(model.integer, dtype=bool)]
        assert np.max(np.abs(integer_values - np.rint(integer_values))) < 1e-4
        assert min(solution - model.lb) >= -1e-4 and min(model.ub - solution) >= -1e-4
        for key, vector in objectives(model).items():
            if key != 'cost':
                value = round(float(-vector @ solution))
                assert row.get(f'{key}_value') in (None, value)
                row[f'final_{key}'] = value
        assert metrics['critical_unserved_kwh'] + 1e-5 >= physical_floor
        row.update(metrics)
        prefix = stem.relative_to(P).as_posix()
        if replay:
            archive = read(prefix + '_replay.csv')
            assert len(trace) == len(archive)
            for field in ('min_energized_voltage', 'max_energized_voltage', 'max_current_ratio',
                          'energized_service_fraction'):
                assert np.allclose([float(r[field]) for r in archive], [r[field] for r in trace], atol=1e-6)
        else:
            if solved is not None:
                assert metrics['AC_voltage_thermal_feasible'], stem
                log_path.write_text(json.dumps(logs, indent=2), encoding='utf-8')
            np.savez_compressed(vector_path, x=solution)
            write(prefix + '_replay.csv', trace)
            write(prefix + '_current.csv', currents)
            write(prefix + '_starts.csv', starts)
    result_path = stem.with_suffix('.result.json')
    if replay:
        saved = json.loads(result_path.read_text(encoding='utf-8'))
        assert saved.keys() == row.keys(), stem
        for key, value in saved.items():
            if isinstance(value, float):
                assert np.isclose(value, row[key], atol=1e-6, rtol=1e-8), (stem, key)
            else:
                assert value == row[key], (stem, key)
    else:
        result_path.write_text(json.dumps(row, indent=2), encoding='utf-8')
    print(json.dumps(row), flush=True)
    return row


def certify_tradeoff(event, seconds, replay):
    settings, _ = settings_for(event, 'long_reserve')
    model = Model(settings, 'coordinated', linear_data(settings), tighten_v=.0003, tighten_i=.012)
    add_island_cut(model, 18, 320)
    care = DIRECTORY / f'{event}_long_reserve_coordinated'
    energy = DIRECTORY / f'{event}_long_reserve_energy_first'
    care_row, energy_row = [json.loads(stem.with_suffix('.result.json').read_text(encoding='utf-8'))
                            for stem in (care, energy)]
    assert care_row['essential_certified'] and care_row['urgent_certified']
    assert care_row['final_essential'] == energy_row['final_essential']
    assert care_row['AC_voltage_thermal_feasible'] and energy_row['AC_voltage_thermal_feasible']
    criteria = objectives(model)
    held = ['essential', 'urgent']
    model.M = vstack([model.M, csr_matrix(np.stack([criteria[k] for k in held]))], format='csc')
    levels = [-care_row[f'final_{k}'] for k in held]
    model.lo.extend(levels)
    model.hi.extend(levels)
    archive = DIRECTORY / f'{event}_long_reserve_certificate'
    if replay:
        logs = json.loads(archive.with_suffix('.json').read_text(encoding='utf-8'))
        vector = np.load(archive.with_suffix('.npz'))['x']
    else:
        vector, logs = solve(model, seconds, archive.with_suffix('.json'), care.with_suffix('.npz'), stages=('energy',), prime=True)
        np.savez_compressed(archive.with_suffix('.npz'), x=vector)
    assert np.max(np.asarray(model.lo) - model.M @ vector) < 1e-4
    assert np.max(model.M @ vector - np.asarray(model.hi)) < 1e-4
    integers = vector[np.asarray(model.integer, dtype=bool)]
    assert np.max(np.abs(integers - np.rint(integers))) < 1e-4
    assert min(vector - model.lb) >= -1e-4 and min(model.ub - vector) >= -1e-4
    assert logs[0]['value'] == round(float(-criteria['energy'] @ vector))
    cap = logs[0]['upper_bound']
    metrics, _, _, _, _, _ = evaluate(settings, model, vector)
    result = dict(event=event, urgent_maximum=care_row['final_urgent'],
                  energy_candidate_kwh=energy_row['final_energy'], upper_at_urgent_max_kwh=cap,
                  treatment_energy_at_urgent_max_kwh=logs[0]['value'],
                  conditional_AC_pass=metrics['AC_voltage_thermal_feasible'],
                  proven=cap is not None and cap < energy_row['final_energy'])
    print(result, flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--event', choices=EVENTS)
    parser.add_argument('--case', choices=CASES)
    parser.add_argument('--policy', choices=('coordinated', 'voltage_only', 'energy_first'))
    parser.add_argument('--seconds', type=int, default=50)
    parser.add_argument('--reserve', action='store_true', help='Size the essential-only fuel stock')
    parser.add_argument('--certificate', action='store_true', help='Bound energy at the certified urgent maximum')
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    DIRECTORY.mkdir(exist_ok=True)
    events = (args.event,) if args.event else EVENTS
    if args.certificate:
        rows = [certify_tradeoff(event, args.seconds, args.replay) for event in events]
        if not args.event:
            write('Severe_Tradeoff_Results.csv', rows)
        return
    if args.reserve:
        names = (args.case,) if args.case else ('mid100', 'long100')
        rows = []
        for event in events:
            for name in names:
                settings, _ = settings_for(event, name)
                stem = DIRECTORY / f'{event}_{name}_reserve_size'
                result = size_reserve(settings, stem, args.replay, 400.)
                result['isolation_h'] = settings['tie_delay_h']
                rows.append(result)
        if not args.event and not args.case:
            write('Severe_Reserve_Results.csv', rows)
        return
    names = (args.case,) if args.case else CASES
    rows = [run(event, name, policy, args.seconds, args.replay) for event in events for name in names
            for policy in ((args.policy,) if args.policy else POLICIES[name])]
    if not args.event and not args.case and not args.policy:
        write('Severe_Stress_Results.csv', [r for r in rows if r['case'] != 'long_no_surge'])
    if not args.event and args.case in (None, 'long_no_surge') and args.policy in (None, 'voltage_only'):
        write('No_Surge_Results.csv', [r for r in rows if r['case'] == 'long_no_surge'])


if __name__ == '__main__':
    main()
