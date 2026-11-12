"""Rule out one more on-time urgent start at a certified energy level."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack

from certify_energy_first import add_urgent_totals, group_starts
from compound_stress import EVENTS, objectives, solve
from severe_stress import Model, add_island_cut, linear_data, settings_for


def check_classes(event, energy, urgent, folder, seconds, replay):
    settings, _ = settings_for(event, 'long_reserve')
    model = Model(settings, 'energy_first', linear_data(settings), tighten_v=.0003, tighten_i=.012)
    add_island_cut(model, 18, 320)
    for index in model.names['plant'].ravel():
        model.lb[index] = model.ub[index] = 1.
    first = model.n
    groups = group_starts(model)
    totals = np.zeros((2, model.n))
    for index, (_, kind, _) in enumerate(groups, start=first):
        totals[kind, index] = 1.
        model.ub[index] = min(model.ub[index], (10, 20)[kind])

    criteria = objectives(model)
    assert all(criteria['energy'][index] == -(26, 3)[kind]
               for index, _, kind, *_ in model.jobs)
    lower, upper = np.asarray(model.lo), np.asarray(model.hi)
    equal = lower == upper
    has_upper = np.isfinite(upper) & ~equal
    has_lower = np.isfinite(lower) & ~equal
    inequalities = vstack([model.M[has_upper], -model.M[has_lower]], format='csc')
    limits = np.r_[upper[has_upper], -lower[has_lower]]
    equalities = vstack([model.M[equal], csr_matrix(totals)], format='csc')
    available_urgent, available_other = settings['A'].sum(axis=(0, 2))
    checks = []
    for high in range(int(urgent), int(available_urgent) + 1):
        other, remainder = divmod(energy - 26 * high, 3)
        if remainder or not 0 <= other <= available_other:
            continue
        relaxation = linprog(criteria['urgent'], A_ub=inequalities, b_ub=limits,
                             A_eq=equalities, b_eq=np.r_[upper[equal], high, other],
                             bounds=np.column_stack([model.lb, model.ub]), method='highs-ds')
        assert relaxation.status in (0, 2), relaxation.message
        check = dict(urgent_patients=high, other_patients=other, relaxation_status=int(relaxation.status))
        if relaxation.success:
            check['relaxed_urgent_upper'] = float(-relaxation.fun)
        checks.append(check)

    add_urgent_totals(model, groups, np.zeros(model.n), criteria)
    totals = np.pad(totals, ((0, 0), (0, model.n - totals.shape[1])))
    model.M = vstack([model.M, csr_matrix(np.vstack([totals, criteria['urgent']]))], format='csc')
    model.lo.extend([0., 0., -np.inf])
    model.hi.extend([0., 0., -urgent])
    for check in checks:
        if check['relaxation_status'] == 2 or check['relaxed_urgent_upper'] < urgent - 1e-6:
            continue
        high, other = check['urgent_patients'], check['other_patients']
        model.lo[-3:-1] = model.hi[-3:-1] = [high, other]
        log = folder / f'class-{high}.json'
        if replay:
            logs = json.loads(log.read_text(encoding='utf-8'))
        else:
            _, logs = solve(model, seconds, log, None, stages=('urgent',), criteria=criteria)
        check['integer_search'] = logs[0]
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('event', choices=EVENTS)
    parser.add_argument('energy', type=int)
    parser.add_argument('urgent', type=int)
    parser.add_argument('folder', type=Path)
    parser.add_argument('--seconds', type=int, default=350)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    args.folder.mkdir(parents=True, exist_ok=True)
    checks = check_classes(args.event, args.energy, args.urgent, args.folder, args.seconds, args.replay)
    (args.folder / 'class_checks.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
    assert all(check.get('integer_search', {}).get('status', 8) == 8 for check in checks), checks
    print(f'No schedule starts at least {args.urgent} urgent patients on time at {args.energy} kWh.')


if __name__ == '__main__':
    main()
