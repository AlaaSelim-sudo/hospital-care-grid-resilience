"""Minimize diesel after fixing essential service and on-time urgent starts."""

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, vstack

from certify_energy_first import add_urgent_totals, group_starts, recover_cohorts
from compound_stress import EVENTS, objectives, solve
from severe_stress import DIRECTORY, Model, P, add_island_cut, evaluate, linear_data, settings_for, write


OUTPUT = P / 'care_diesel_data'
CASES = ('long_no_surge', 'long_reserve')


def diesel_objective(model):
    objective = np.zeros(model.n)
    objective[model.names['generator']] = .30
    objective[model.names['generator_on']] = 1.5
    return objective


def validate(model, solution):
    result = model.M @ solution
    assert np.all(result >= np.asarray(model.lo) - 1e-5)
    assert np.all(result <= np.asarray(model.hi) + 1e-5)
    assert np.all(solution >= np.asarray(model.lb) - 1e-5)
    assert np.all(solution <= np.asarray(model.ub) + 1e-5)
    integers = solution[np.asarray(model.integer, dtype=bool)]
    assert np.max(np.abs(integers - np.rint(integers))) < 1e-5


def run(event, case, seconds):
    settings, _ = settings_for(event, case)
    model = Model(settings, 'coordinated', linear_data(settings), tighten_v=.0003, tighten_i=.012)
    add_island_cut(model, 18, 320)
    essential_slots = model.names['plant'].size
    stem = OUTPUT / f'{event}_{case}'
    witness = stem.with_suffix('.npz')
    with tempfile.TemporaryDirectory(dir=OUTPUT) as folder:
        work = Path(folder)
        if witness.exists():
            initial = np.load(witness)['x']
            previous = json.loads(stem.with_suffix('.json').read_text(encoding='utf-8'))
            priorities = previous[:2]
        elif case == 'long_reserve':
            saved = DIRECTORY / f'{event}_{case}_coordinated'
            initial = np.load(saved.with_suffix('.npz'))['x']
            priorities = json.loads(saved.with_suffix('.json').read_text(encoding='utf-8'))[:2]
        else:
            initial, priorities = solve(model, seconds, work / 'care.json', None,
                                        stages=('essential', 'urgent'))
        assert initial is not None
        assert [r['stage'] for r in priorities] == ['essential', 'urgent']
        assert all(r['certified'] for r in priorities), priorities
        assert priorities[0]['value'] == essential_slots
        urgent = priorities[1]['value']
        for idx in model.names['plant'].ravel():
            model.lb[idx] = model.ub[idx] = 1.
        physical_model = Model(settings, 'coordinated', linear_data(settings), tighten_v=.0003, tighten_i=.012)
        add_island_cut(physical_model, 18, 320)
        validate(physical_model, initial)
        old_size = model.n
        groups = group_starts(model)
        start_totals = [sum(initial[i] for i in members) for members in groups.values()]
        warm = np.r_[initial, start_totals]
        criteria = objectives(model)
        criteria['diesel'] = diesel_objective(model)
        for offset, (_, kind, _) in enumerate(groups, start=old_size):
            model.ub[offset] = min(model.ub[offset], (10, 20)[kind])
        warm = add_urgent_totals(model, groups, warm, criteria)
        model.M = vstack([model.M, csr_matrix(criteria['urgent'])], format='csc')
        model.lo.append(-urgent)
        model.hi.append(-urgent)
        seed = work / 'seed.npz'
        np.savez_compressed(seed, x=warm)
        solution, logs = solve(model, seconds, work / 'diesel.json', seed,
                               stages=('diesel',), prime=True, criteria=criteria)
        assert solution is not None
        validate(model, solution)
        totals = np.rint(solution[old_size:old_size + len(groups)])
        original = solution[:old_size].copy()
        recover_cohorts(physical_model.jobs, original, dict(zip(groups, totals)), settings['A'])
        validate(physical_model, original)
        metrics, trace, _, currents, starts, dispatch = evaluate(settings, physical_model, original)
        assert metrics['high_timely'] == urgent and metrics['critical_unserved_kwh'] == 0
        assert metrics['AC_voltage_thermal_feasible'], metrics
        diesel = float(diesel_objective(physical_model) @ original)
        bound = logs[0]['dual_bound']
        if witness.exists():
            old_bound = previous[2]['lower_bound']
            if old_bound is not None:
                bound = old_bound if bound is None else max(old_bound, bound)
        assert np.isclose(diesel, metrics['fuel_litres'])
        assert bound is None or bound <= diesel + 1e-5
        gap = None if bound is None else max(0., diesel - bound)
        logs[0].update(value=diesel, lower_bound=bound, absolute_gap_l=gap,
                       certified=gap is not None and gap <= 1e-5)
        np.savez_compressed(witness, x=original)
        stem.with_suffix('.json').write_text(json.dumps([*priorities, *logs], indent=2), encoding='utf-8')
        result = dict(event=event, case=case, essential_site_hours=essential_slots,
                      diesel_lower_bound_l=bound, diesel_gap_l=gap, diesel_certified=logs[0]['certified'],
                      **metrics)
        stem.with_suffix('.result.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
        prefix = stem.relative_to(P).as_posix()
        write(prefix + '_replay.csv', trace)
        write(prefix + '_current.csv', currents)
        write(prefix + '_starts.csv', starts)
        rows = []
        for h in range(3):
            for t in range(dispatch['generator'].shape[1]):
                generation = float(dispatch['generator'][h, t])
                rows.append(dict(site=f'H{h+1}', hour=t, generator_kw=generation,
                                 diesel_l=.30*generation + 1.5*float(dispatch['generator_on'][h, t]),
                                 battery_discharge_kw=float(dispatch['discharge'][h, t])))
        write(prefix + '_resources.csv', rows)
        print(json.dumps(result), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('event', choices=EVENTS)
    parser.add_argument('case', choices=CASES)
    parser.add_argument('--seconds', type=int, default=180)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    run(args.event, args.case, args.seconds)


if __name__ == '__main__':
    main()
