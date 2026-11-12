"""Matched care-first and energy-first tests under tighter transformer ratings."""

import argparse
import json

from run_fault_study import Model, P, evaluate, linear_data, np, scenario, write


def compare(event, rating, seconds):
    settings = scenario(event, rating=rating)
    linear = linear_data(settings)
    directory = P / 'objective_comparison'
    directory.mkdir(exist_ok=True)
    rows = []
    for policy in ('energy_first', 'coordinated'):
        print(f'Solving {event}, {rating} kVA, {policy}', flush=True)
        # The two policies share the same linearization and guard bands.
        model = Model(settings, policy, linear, tighten_v=.0003, tighten_i=.012)
        solution, stages = model.optimize(seconds)
        row = dict(event=event, rating_kva=rating, policy=policy, solve_limit_s=seconds)
        stem = directory / f'{event}_{rating}_{policy}'
        stem.with_suffix('.json').write_text(json.dumps(stages, indent=2))
        row['primary_gap'] = stages[1].get('gap') if len(stages) > 1 else None
        row['essential_gap'] = stages[0].get('gap')
        row['clinical_stage_found'] = any(
            stage['stage'] in ('clinical_tie_break', 'care') and stage['objective'] is not None
            for stage in stages
        )
        if solution is not None:
            metrics, trace, _, _, starts, _ = evaluate(settings, model, solution)
            row.update(metrics)
            np.savez_compressed(stem.with_suffix('.npz'), x=solution)
            write(stem.relative_to(P).as_posix() + '_starts.csv', starts)
            write(stem.relative_to(P).as_posix() + '_replay.csv', trace)
        rows.append(row)
        write(f'objective_comparison/{event}_{rating}.csv', rows)
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('event', choices=['dec2023', 'apr2024'])
    parser.add_argument('rating', type=int, choices=[500, 550, 600])
    parser.add_argument('--seconds', type=int, default=60, help='MILP limit per objective stage')
    args = parser.parse_args()
    compare(args.event, args.rating, args.seconds)
