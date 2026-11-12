"""Recheck the archived capacity-test schedules with OpenDSS."""

import math

from run_fault_study import Model, P, evaluate, linear_data, np, read, scenario


def check_schedule(path, settings, policy, recorded, fixed_care=None):
    model = Model(settings, policy, linear_data(settings), tighten_v=.0003, tighten_i=.012)
    with np.load(path) as archive:
        solution = archive['x']
    integer = np.asarray(model.integer, dtype=bool)
    assert np.max(abs(solution[integer] - np.rint(solution[integer]))) < 1e-5
    assert np.all(solution >= np.asarray(model.lb) - 1e-5)
    assert np.all(solution <= np.asarray(model.ub) + 1e-5)
    if fixed_care is not None:
        indices = [*model.names['plant'].ravel(), *(job[0] for job in model.jobs)]
        assert np.array_equal(np.rint(solution[indices]), np.rint(fixed_care[indices]))
    actual, _, _, loading, _, _ = evaluate(settings, model, solution)
    assert actual['linear_residual'] < 1e-5
    for key in ('timely', 'high_timely', 'critical_unserved_kwh', 'voltage_violation_h',
                'overload_h', 'max_current_ratio', 'fuel_litres'):
        assert math.isclose(actual[key], float(recorded[key]), rel_tol=1e-6, abs_tol=1e-6), (path, key)
    assert str(actual['AC_voltage_thermal_feasible']) == recorded['AC_voltage_thermal_feasible']
    if recorded.get('tie_peak_current_a'):
        peak = max(item['max_current_ratio'] for item in loading
                   if item['element'].lower() == 'line.restorationtie')
        assert math.isclose(peak * settings['tie_rating_A'], float(recorded['tie_peak_current_a']),
                            rel_tol=1e-6)


def main():
    checked = 0
    for event in ('dec2023', 'apr2024'):
        for rating in (500, 550):
            for row in read(f'objective_comparison/{event}_{rating}.csv'):
                if not row.get('AC_voltage_thermal_feasible'):
                    continue
                policy = row['policy']
                check_schedule(P / 'objective_comparison' / f'{event}_{rating}_{policy}.npz',
                               scenario(event, rating=rating), policy, row)
                checked += 1
    for row in read('Tie_Redispatch_Results.csv'):
        if not row.get('AC_voltage_thermal_feasible'):
            continue
        event, rating = row['case'], int(row['tie_rating_a'])
        with np.load(P / 'solutions' / f'v4_{event}_r12_d6_x600_a500_coordinated.npz') as archive:
            original = archive['x']
        check_schedule(P / 'tie_redispatch' / f'{event}_{rating}.npz',
                       scenario(event, tie_rating=rating), 'coordinated', row, original)
        checked += 1
    print(f'Validated {checked} archived capacity-test schedules against AC replay.')


if __name__ == '__main__':
    main()
