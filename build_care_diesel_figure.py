"""Create the matched diesel, urgent care, and feeder-loading TikZ panels."""

import json

import numpy as np

from build_load_figure import annotation, line
from minimize_care_diesel import OUTPUT
from run_fault_study import P, read, write
from severe_stress import settings_for


EVENTS = (('dec2023', 'December 2023'), ('apr2024', 'April 2024'))
CASES = (('long_no_surge', 'teal,dashed', 'No surge'),
         ('long_reserve', 'navy,solid', 'Surge'))


def load_result(event, case, deadline_curve):
    stem = OUTPUT / f'{event}_{case}'
    result = json.loads(stem.with_suffix('.result.json').read_text(encoding='utf-8'))
    stem_name = stem.relative_to(P).as_posix()
    settings, _ = settings_for(event, case)
    fault = settings['fault_hour']
    assert (fault, settings['tie_delay_h'], settings['repair_h']) == (12, 18, 27)
    assert result['essential_site_hours'] == 216 and result['AC_voltage_thermal_feasible']
    logs = json.loads(stem.with_suffix('.json').read_text(encoding='utf-8'))
    assert [stage['stage'] for stage in logs] == ['essential', 'urgent', 'diesel']
    assert all(stage['certified'] for stage in logs[:2])
    assert logs[1]['value'] == result['high_timely']
    fuel_rows = read(stem_name + '_resources.csv')
    by_hour = np.zeros(len(settings['stages']))
    for row in fuel_rows:
        by_hour[int(row['hour'])] += float(row['diesel_l'])
    assert np.isclose(by_hour.sum(), result['fuel_litres'])
    fuel = list(enumerate(np.r_[0., np.cumsum(by_hour)], start=-fault))
    replay = read(stem_name + '_replay.csv')
    current = [(float(row['hour']) - fault, 100 * float(row['max_current_ratio'])) for row in replay]
    assert [hour for hour, _ in current] == [i / 4 - fault for i in range(288)]
    assert np.isclose(max(value for _, value in current), 100 * result['max_current_ratio'])
    clinical = deadline_curve(read(stem_name + '_starts.csv'), settings,
                              {**result, 'AC_voltage_thermal_feasible': 'True'})
    return result, {'fuel': fuel, 'missed': [*clinical, (60, clinical[-1][1])],
                    'current': [*current, (60, current[-1][1])]}


def diesel_figure(deadline_curve_fn):
    # The deadline definition is shared with the existing severe-case figures.
    data = {(event, case): load_result(event, case, deadline_curve_fn)
            for event, _ in EVENTS for case, _, _ in CASES}
    summary = [record for record, _ in data.values()]
    write('Care_Diesel_Results.csv', summary)

    body = [r'''\providecommand{\DieselPanelHeight}{1.68cm}
\providecommand{\DieselPanelGap}{.50cm}
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 3,horizontal sep=.86cm,vertical sep=\DieselPanelGap,
 xlabels at=edge bottom,xticklabels at=edge bottom},
 width=.415\textwidth,height=\DieselPanelHeight,scale only axis,xmin=-12,xmax=60,
 xtick={-12,0,12,18,27,42,60},xlabel={Hours from fault},
 tick label style={font=\scriptsize},label style={font=\scriptsize},
 title style={font=\footnotesize\bfseries},ymajorgrids,grid style={gray!20},
 axis line style={black!65}]
''']
    panels = [('fuel', 'Diesel (L)', 0, 700, '0,300,600'),
              ('missed', 'Misses', 0, 50, '0,20,40'),
              ('current', 'Peak I (\\%)', 76, 105, '80,90,100')]
    for field, axis_label, low, high, ticks in panels:
        for event, title in EVENTS:
            options = [f'ymin={low},ymax={high},ytick={{{ticks}}}']
            if field == 'fuel':
                options.append(f'title={{{title}}}')
            if event == EVENTS[0][0]:
                options.append(f'ylabel={{{axis_label}}}')
            body.append(r'\nextgroupplot[' + ','.join(options) + ']')
            body.append(rf'\fill[navy!6] (axis cs:0,{low}) rectangle (axis cs:18,{high});')
            for hour in (0, 18, 27):
                body.append(rf'\draw[black!52,densely dashed] (axis cs:{hour},{low})--(axis cs:{hour},{high});')
            if field == 'current':
                body.append(r'\draw[black!78,densely dotted] (axis cs:-12,100)--(axis cs:60,100);')
            for case, style, _ in CASES:
                _, traces = data[event, case]
                body.append(line(style, traces[field], stepwise=field != 'fuel'))
            base, surge = [data[event, case][0] for case, _, _ in CASES]
            if field == 'fuel':
                label = rf'{base["fuel_litres"]:.1f} $\to$ {surge["fuel_litres"]:.1f} L'
                body.append(annotation(59, 70, label, 'anchor=east,fill=white,inner sep=1pt'))
            elif field == 'missed':
                totals = [f'{r["high_arrivals"]-r["high_timely"]}/{r["high_arrivals"]}' for r in (base, surge)]
                body.append(annotation(59, 7, rf'{totals[0]} $\to$ {totals[1]}',
                                       'anchor=east,fill=white,inner sep=1pt'))
            else:
                label = rf'{100*base["max_current_ratio"]:.1f}\% $\to$ {100*surge["max_current_ratio"]:.1f}\%'
                body.append(annotation(59, 77.5, label, 'anchor=east,fill=white,inner sep=1pt'))
    body.append(r'\end{groupplot}')
    keys = [rf'\tikz[baseline=-.45ex]\draw[{style},line width=.8pt] (0,0)--(.35,0);\,{label}'
            for _, style, label in CASES]
    body.append(r'\node[anchor=south,font=\scriptsize] at ([xshift=.43\textwidth,yshift=.46cm]group c1r1.north west) {' +
                r'\quad '.join(keys) + '};')
    body.append(r'\end{tikzpicture}')
    (P / 'figures' / 'care_diesel_comparison.tikz').write_text('\n'.join(body), encoding='utf-8')
