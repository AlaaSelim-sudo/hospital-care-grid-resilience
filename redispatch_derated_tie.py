"""Test whether the main care schedule can be retained after tie derating."""

import json

from run_fault_study import Bounds, LinearConstraint, Model, P, evaluate, linear_data, milp, np, scenario, write


def main():
    rows = []
    directory = P / 'tie_redispatch'
    directory.mkdir(exist_ok=True)
    for event in ('dec2023', 'apr2024'):
        saved = np.load(P / 'solutions' / f'v4_{event}_r12_d6_x600_a500_coordinated.npz')
        original = saved['x']
        for rating in (500, 485, 475, 450):
            settings = scenario(event, tie_rating=rating)
            model = Model(settings, 'coordinated', linear_data(settings),
                          tighten_v=.0003, tighten_i=.012)
            fixed = np.concatenate((model.names['plant'].ravel(),
                                    np.array([job[0] for job in model.jobs], dtype=int)))
            lower, upper = np.array(model.lb), np.array(model.ub)
            lower[fixed] = upper[fixed] = np.rint(original[fixed])
            objective = np.zeros(model.n)
            objective[model.names['generator']] = .30
            objective[model.names['generator_on']] = 1.5
            result = milp(objective, integrality=model.integer,
                          bounds=Bounds(lower, upper),
                          constraints=LinearConstraint(model.M, model.lo, model.hi),
                          options={'time_limit': 45, 'mip_rel_gap': 1e-5})
            row = dict(case=event, tie_rating_a=rating,
                       linear_solver_status=int(result.status), solver_message=result.message)
            if result.x is not None:
                metrics, trace, _, loading, _, _ = evaluate(settings, model, result.x)
                tie_peak = max(item['max_current_ratio'] for item in loading
                               if item['element'].lower() == 'line.restorationtie')
                row.update(metrics, tie_peak_current_a=rating * tie_peak)
                stem = directory / f'{event}_{rating}'
                np.savez_compressed(stem.with_suffix('.npz'), x=result.x)
                write(stem.relative_to(P).as_posix() + '_replay.csv', trace)
            rows.append(row)
            write('Tie_Redispatch_Results.csv', rows)
            print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
