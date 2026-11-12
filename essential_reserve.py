"""Size a common pre-event fuel stock for essential-only service in the extended case."""

import argparse
import json

from scipy.sparse import coo_matrix, hstack, vstack

from compound_stress import DIRECTORY, EVENTS, design
from run_fault_study import Bounds, LinearConstraint, Model, evaluate, linear_data, milp, np, write


def run(event, replay):
    settings, _, _ = design(event, 'extended')
    archive = DIRECTORY / f'{event}_essential_reserve'
    return size_reserve(settings, archive, replay, 300.)


def size_reserve(settings, archive, replay, maximum_stock):
    settings = settings.copy()
    settings['fuel_limit_l'] = (np.inf, np.inf, 160.)
    model = Model(settings, 'coordinated', linear_data(settings), tighten_v=.0003, tighten_i=.012)
    lower, upper = np.array(model.lb + [0.]), np.array(model.ub + [maximum_stock])
    lower[model.names['plant']] = 1
    upper[[job[0] for job in model.jobs]] = 0
    # One tank capacity bounds each islanded site's own fuel use; it is not a pooled budget.
    fuel_rows = np.zeros((2, model.n + 1))
    for h, row in enumerate(fuel_rows):
        row[model.names['generator'][h]] = .30
        row[model.names['generator_on'][h]] = 1.5
        row[-1] = -1
    matrix = vstack([hstack([model.M, coo_matrix((model.M.shape[0], 1))]), coo_matrix(fuel_rows)]).tocsc()
    lo, hi = np.array(model.lo + [-np.inf] * 2), np.array(model.hi + [0.] * 2)
    objective = np.zeros(model.n + 1)
    objective[-1] = 1
    if replay:
        x = np.load(archive.with_suffix('.npz'))['x']
        log = json.loads(archive.with_suffix('.json').read_text(encoding='utf-8'))
    else:
        result = milp(objective, integrality=model.integer + [0], bounds=Bounds(lower, upper),
                      constraints=LinearConstraint(matrix, lo, hi), options={'time_limit': 75, 'mip_rel_gap': 1e-5})
        assert result.x is not None, result.message
        x = result.x
        log = model.log(result, 'common_initial_fuel_l')
        np.savez_compressed(archive.with_suffix('.npz'), x=x)
        archive.with_suffix('.json').write_text(json.dumps(log, indent=2), encoding='utf-8')
    assert np.max(lo - matrix @ x) < 1e-5 and np.max(matrix @ x - hi) < 1e-5
    assert np.max(lower - x) < 1e-5 and np.max(x - upper) < 1e-5
    integer_values = x[:-1][np.asarray(model.integer, dtype=bool)]
    assert np.max(np.abs(integer_values - np.rint(integer_values))) < 1e-5
    metrics, _, _, _, _, _ = evaluate(settings, model, x[:-1])
    assert metrics['completed'] == 0 and metrics['critical_unserved_kwh'] < 1e-5
    row = dict(event=settings['case'], capacity_per_island_l=x[-1], lower_bound_l=log['dual_bound'],
               solver_gap=log['gap'], AC_pass=metrics['AC_voltage_thermal_feasible'])
    print(row, flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    DIRECTORY.mkdir(exist_ok=True)
    write('Essential_Reserve_Results.csv', [run(event, args.replay) for event in EVENTS])


if __name__ == '__main__':
    main()
