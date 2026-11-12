"""Certify the first three priorities for the matched 18-hour energy baseline."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csr_matrix, hstack, vstack

from compound_stress import EVENTS, objectives, solve
from severe_stress import (DIRECTORY, Model, add_island_cut, certify_tradeoff, evaluate,
                           linear_data, read, run, settings_for, write)


def group_starts(model):
    """Use integer totals by site, class and start time with cohort assignment flows."""
    groups = defaultdict(list)
    for index, site, kind, _, start, _ in model.jobs:
        groups[site, kind, start].append(index)

    old_size = model.n
    old_rows = model.M.shape[0]
    count = len(groups)
    row_indices, col_indices, coefficients = [], [], []
    for row, members in enumerate(groups.values()):
        row_indices.extend([row] * (len(members) + 1))
        col_indices.extend([*members, old_size + row])
        coefficients.extend([1.] * len(members) + [-1.])
        for member in members:
            model.integer[member] = 0
    # Use the integer totals directly in rows where every assignment in a
    # group has the same coefficient. Cohort-capacity rows keep their flows.
    group_for = {member: (old_size + row, len(members))
                 for row, members in enumerate(groups.values()) for member in members}
    original = model.M.tocsr()
    main_rows, main_cols, main_vals = [], [], []
    for row in range(old_rows):
        start, stop = original.indptr[row:row + 2]
        pending = defaultdict(list)
        for column, value in zip(original.indices[start:stop], original.data[start:stop]):
            group = group_for.get(column)
            if group is None:
                main_rows.append(row)
                main_cols.append(column)
                main_vals.append(value)
            else:
                pending[group].append((column, value))
        for (column, size), entries in pending.items():
            coefficient = entries[0][1]
            terms = [(column, coefficient)] if (
                len(entries) == size and all(value == coefficient for _, value in entries)
            ) else entries
            for target, value in terms:
                main_rows.append(row)
                main_cols.append(target)
                main_vals.append(value)
    grouped = coo_matrix((main_vals, (main_rows, main_cols)), shape=(old_rows, old_size + count))
    linking = coo_matrix((coefficients, (row_indices, col_indices)),
                         shape=(count, old_size + count))
    model.M = vstack([grouped, linking], format='csc')
    model.lo.extend([0.] * count)
    model.hi.extend([0.] * count)
    model.lb.extend([0.] * count)
    model.ub.extend(sum(model.ub[i] for i in members) for members in groups.values())
    model.integer.extend([1] * count)
    model.n += count
    return groups


def recover_cohorts(jobs, vector, group_totals, arrival_counts):
    """Recover integral cohort assignments using a bipartite flow LP."""
    by_site_class = defaultdict(list)
    for index, site, kind, arrival, start, _ in jobs:
        by_site_class[site, kind].append((index, arrival, start))
    for (site, kind), edges in by_site_class.items():
        arrivals = {hour: i for i, hour in enumerate(sorted({edge[1] for edge in edges}))}
        starts = {hour: i for i, hour in enumerate(sorted({edge[2] for edge in edges}))}
        columns = np.arange(len(edges))
        ones = np.ones(len(edges))
        inbound = csr_matrix((ones, ([arrivals[r] for _, r, _ in edges], columns)),
                             shape=(len(arrivals), len(edges)))
        outbound = csr_matrix((ones, ([starts[s] for _, _, s in edges], columns)),
                              shape=(len(starts), len(edges)))
        reward = np.asarray([-(kind == 0 and s - r <= 2) for _, r, s in edges], dtype=float)
        assignment = linprog(reward, A_ub=inbound, b_ub=[arrival_counts[site, kind, r] for r in arrivals],
                             A_eq=outbound, b_eq=[group_totals[site, kind, s] for s in starts],
                             method='highs-ds')
        assert assignment.success, assignment.message
        counts = np.rint(assignment.x)
        assert np.max(np.abs(counts - assignment.x)) < 1e-5
        vector[[i for i, _, _ in edges]] = counts


def add_urgent_totals(model, groups, initial, criteria):
    """Let the solver branch on integer counts of on-time starts."""
    timely = [[i for i in members if criteria['urgent'][i] == -1]
              for members in groups.values()]
    timely = [members for members in timely if members]
    first = model.n
    count = len(timely)
    rr, cc, values = [], [], []
    for row, members in enumerate(timely):
        rr.extend([row] * (len(members) + 1))
        cc.extend([*members, first + row])
        values.extend([1.] * len(members) + [-1.])
        criteria['urgent'][members] = 0.
    linking = coo_matrix((values, (rr, cc)), shape=(count, first + count))
    model.M = vstack([hstack([model.M, csr_matrix((model.M.shape[0], count))]),
                      linking], format='csc')
    model.lo.extend([0.] * count)
    model.hi.extend([0.] * count)
    model.lb.extend([0.] * count)
    model.ub.extend(sum(model.ub[i] for i in members) for members in timely)
    model.integer.extend([1] * count)
    model.n += count
    for name, objective in criteria.items():
        criteria[name] = np.r_[objective, np.full(count, -1. if name == 'urgent' else 0.)]
    return np.r_[initial, [sum(initial[i] for i in members) for members in timely]]


def certify(event, folder, seconds, stage, seed=None, energy_target=None, urgent_groups=False,
            min_urgent=None):
    settings, _ = settings_for(event, 'long_reserve')
    model = Model(settings, 'energy_first', linear_data(settings), tighten_v=.0003, tighten_i=.012)
    add_island_cut(model, 18, 320)
    # A saved feasible solution serves all 216 site-hours, the physical upper bound.
    for index in model.names['plant'].ravel():
        model.lb[index] = model.ub[index] = 1.
    source = seed or DIRECTORY / f'{event}_long_reserve_energy_first.npz'
    vector = np.load(source)['x']
    old_size = model.n
    groups = group_starts(model)
    initial = np.r_[vector, [sum(vector[i] for i in members) for members in groups.values()]]
    criteria = objectives(model)
    energy = criteria['energy']
    for offset, members in enumerate(groups.values(), start=old_size):
        assert np.all(energy[members] == energy[members[0]])
        energy[offset] = energy[members[0]]
        energy[members] = 0.
    for offset, (_, kind, _) in enumerate(groups, start=old_size):
        model.ub[offset] = min(model.ub[offset], (10, 20)[kind])
    if urgent_groups:
        initial = add_urgent_totals(model, groups, initial, criteria)
        energy = criteria['energy']
    if stage == 'urgent':
        assert round(-energy @ initial) == energy_target
        model.M = vstack([model.M, csr_matrix(energy)], format='csc')
        model.lo.append(-energy_target)
        model.hi.append(-energy_target)
    if min_urgent is not None:
        assert round(-criteria['urgent'] @ initial) >= min_urgent
        model.M = vstack([model.M, csr_matrix(criteria['urgent'])], format='csc')
        model.lo.append(-np.inf)
        model.hi.append(-min_urgent)
    criterion = criteria[stage]
    model.M = vstack([model.M, csr_matrix(criterion)], format='csc')
    model.lo.append(-np.inf)
    model.hi.append(round(criterion @ initial))
    stem = folder / f'{event}_{stage}'
    warm = stem.with_suffix('.seed.npz')
    np.savez_compressed(warm, x=initial)
    solution, logs = solve(model, seconds, stem.with_suffix('.json'), warm,
                           stages=(stage,), prime=True, criteria=criteria)
    assert solution is not None
    totals = solution[old_size:old_size + len(groups)]
    group_totals = dict(zip(groups, np.rint(totals)))
    assert np.max(np.abs(totals - np.rint(totals))) < 1e-5
    recover_cohorts(model.jobs, solution, group_totals, settings['A'])
    value = round(-objectives(model)[stage] @ solution)
    log = logs[0]
    assert log['value'] <= value <= log['upper_bound']
    log.update(value=value, objective=-value, certified=value == log['upper_bound'],
               gap=abs(-value - log['dual_bound']) / max(abs(value), 1))
    if stage == 'urgent':
        log['fixed_energy_kwh'] = energy_target
    if min_urgent is not None:
        log['min_urgent'] = min_urgent
    stem.with_suffix('.json').write_text(json.dumps(logs, indent=2), encoding='utf-8')
    np.savez_compressed(stem.with_suffix('.npz'), x=solution[:old_size])
    return log


def publish(folder):
    """Replace the archived baseline with AC-checked schedules and recorded bounds."""
    cases = read('Severe_Stress_Results.csv')
    completed = {}
    for event in EVENTS:
        energy, urgent = [json.loads((folder / f'{event}_{stage}.json').read_text())[0]
                          for stage in ('energy', 'urgent')]
        assert energy['certified'] and urgent['certified'], (event, energy, urgent)
        assert urgent['fixed_energy_kwh'] == energy['value']
        original = DIRECTORY / f'{event}_long_reserve_energy_first.json'
        essential = json.loads(original.read_text())[0]
        assert essential['certified'] and essential['value'] == 216
        vector = np.load(folder / f'{event}_urgent.npz')['x']
        settings, _ = settings_for(event, 'long_reserve')
        model = Model(settings, 'energy_first', linear_data(settings), tighten_v=.0003, tighten_i=.012)
        add_island_cut(model, 18, 320)
        for log in (essential, energy, urgent):
            assert round(-objectives(model)[log['stage']] @ vector) == log['value']
        metrics, *_ = evaluate(settings, model, vector)
        assert metrics['AC_voltage_thermal_feasible'] and metrics['linear_residual'] < 1e-4
        completed[event, 'energy_first'] = (vector, [essential, energy, urgent])
        # Use the archived secondary-energy solve at the proven urgent maximum.
        care = DIRECTORY / f'{event}_long_reserve_coordinated'
        care_logs = json.loads(care.with_suffix('.json').read_text())[:2]
        tie = DIRECTORY / f'{event}_long_reserve_certificate'
        care_logs += json.loads(tie.with_suffix('.json').read_text())
        care_vector = np.load(tie.with_suffix('.npz'))['x']
        for log in care_logs:
            assert round(-objectives(model)[log['stage']] @ care_vector) == log['value']
        metrics, *_ = evaluate(settings, model, care_vector)
        assert metrics['AC_voltage_thermal_feasible'] and metrics['linear_residual'] < 1e-4
        completed[event, 'coordinated'] = (care_vector, care_logs)
    updated = {(event, policy): run(event, 'long_reserve', policy, 0, False, solved=solution)
               for (event, policy), solution in completed.items()}
    write('Severe_Stress_Results.csv', [
        updated.get((row['event'], row['policy']), row) if row['case'] == 'long_reserve' else row
        for row in cases])
    write('Severe_Tradeoff_Results.csv', [certify_tradeoff(event, 0, True) for event in EVENTS])


def urgent_bound_from_energy(energy, inverse, urgent_value):
    """Rule out an urgent count when its energy bound is below the energy optimum."""
    assert energy['certified'] and inverse['upper_bound'] < energy['value']
    upper = inverse['min_urgent'] - 1
    assert urgent_value <= upper
    return dict(stage='urgent', status=inverse['status'], solver=inverse['solver'],
                message='Urgent bound from the constrained energy solve',
                value=urgent_value, objective=-urgent_value, upper_bound=upper,
                dual_bound=-upper, gap=(upper - urgent_value) / max(urgent_value, 1),
                certified=urgent_value == upper, fixed_energy_kwh=energy['value'],
                inverse_energy_proof=inverse)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('event', choices=EVENTS)
    parser.add_argument('folder', type=Path)
    parser.add_argument('stage', choices=('energy', 'urgent'))
    parser.add_argument('--seconds', type=int, default=300)
    parser.add_argument('--seed', type=Path)
    parser.add_argument('--energy-target', type=int)
    parser.add_argument('--urgent-groups', action='store_true')
    parser.add_argument('--min-urgent', type=int,
                        help='Bound energy subject to this many on-time urgent starts')
    args = parser.parse_args()
    if args.stage == 'urgent' and (args.energy_target is None or args.seed is None):
        parser.error('the urgent stage needs the energy optimum and its solution')
    if args.min_urgent is not None and args.stage != 'energy':
        parser.error('--min-urgent applies to the energy stage')
    args.folder.mkdir(parents=True, exist_ok=True)
    certify(args.event, args.folder, args.seconds, args.stage, args.seed, args.energy_target,
            args.urgent_groups, args.min_urgent)


if __name__ == '__main__':
    main()
