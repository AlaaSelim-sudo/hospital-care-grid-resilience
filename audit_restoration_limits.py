"""Attribute reported AC violations and screen saved schedules at lower tie ratings."""

from collections import defaultdict
from run_fault_study import DT, read, write


def main():
    groups = defaultdict(list)
    for row in read('Fault_Current_Loading.csv'):
        groups[row['case'], row['policy']].append(row)

    attribution = []
    tie_screen = []
    for (case, policy), rows in sorted(groups.items()):
        violating_hours = set()
        counts = defaultdict(set)
        tie_currents = {}
        for row in rows:
            hour = float(row['hour'])
            ratio = float(row['max_current_ratio'])
            element = row['element']
            if ratio > 1 + 1e-7:
                violating_hours.add(hour)
                counts[row['stage'], element].add(hour)
            if element.lower() == 'line.restorationtie':
                tie_currents[hour] = 500 * ratio

        reported = next(
            item for item in read('Fault_Results.csv')
            if item['case'] == case and item['policy'] == policy
        )
        if abs(len(violating_hours) * DT - float(reported['overload_h'])) > 1e-8:
            raise ValueError(f'Overload duration does not match the saved AC result: {case}, {policy}')

        for (stage, element), hours in sorted(counts.items()):
            attribution.append(dict(case=case, policy=policy, stage=stage,
                                    element=element, overloaded_h=len(hours) * DT))
        if tie_currents:
            for rating in (500, 475, 450):
                tie_screen.append(dict(
                    case=case, policy=policy, screened_rating_a=rating,
                    peak_current_a=max(tie_currents.values()),
                    backfeed_h=len(tie_currents) * DT,
                    overloaded_h=DT * sum(current > rating + 1e-7
                                          for current in tie_currents.values()),
                ))

    write('Overcurrent_Attribution.csv', attribution)
    write('Tie_Rating_Screen.csv', tie_screen)


if __name__ == '__main__':
    main()
