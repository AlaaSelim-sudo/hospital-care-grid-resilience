"""Verify the archived quarter-hour voltage results used in Fig. 4."""
from pathlib import Path
import csv

HERE = Path(__file__).resolve().parent
CASES = {
    ('dec2023', 'local', 'no'): (18, 9),
    ('dec2023', 'local', 'yes'): (16.75, 8.25),
    ('apr2024', 'local', 'no'): (17.25, 5),
    ('apr2024', 'local', 'yes'): (18.25, 6.75),
}

def read_trace(storm, policy, surge):
    case = 'long_reserve' if surge == 'yes' else 'long_no_surge'
    dispatch = 'autonomous' if policy == 'local' else 'coordinated'
    with (HERE / f'{storm}_{case}_{dispatch}.csv').open(newline='', encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    assert [float(row['hour']) - 12 for row in rows] == [i / 4 - 12 for i in range(288)]
    return {float(row['hour']) - 12: row for row in rows}

def main():
    for storm in ('dec2023', 'apr2024'):
        for policy in ('local', 'joint'):
            for surge in ('no', 'yes'):
                trace = read_trace(storm, policy, surge)
                low = [t for t, row in trace.items() if float(row['min_energized_voltage']) < .95]
                assert all(float(row['max_energized_voltage']) <= 1.05 for row in trace.values())
                assert all(t < 0 or t >= 18 for t in low)
                post_fault = sum(t >= 0 for t in low) / 4
                backfeed = sum(18 <= t < 27 for t in low) / 4
                assert (post_fault, backfeed) == CASES.get((storm, policy, surge), (0, 0))
                print(f'{storm} {policy:5s} surge={surge:3s}  post-fault={post_fault:5.2f} h  backfeed={backfeed:4.2f} h')
    for surge, expected_v, expected_kw in (('no', .9513, 152), ('yes', .9492, 262)):
        sample = read_trace('apr2024', 'local', surge)[22.75]
        assert round(float(sample['min_energized_voltage']), 4) == expected_v
        assert round(float(sample['hospital_grid_kw'])) == expected_kw
    print('April hour 22.75: 0.9513 to 0.9492 p.u.; hospital imports 152 to 262 kW.')

if __name__ == '__main__':
    main()
