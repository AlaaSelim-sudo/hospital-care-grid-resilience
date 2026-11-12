"""Reconstruct and plot the matched 18-hour surge comparison."""

from collections import defaultdict
from math import ceil

import numpy as np

from run_fault_study import Model, P, evaluate, linear_data, read
from severe_stress import DIRECTORY, settings_for


EVENTS = (('dec2023', 'December 2023'), ('apr2024', 'April 2024'))
SURGE = 'orange!64!black,solid'
BASE = 'teal,dashed'
CARE = 'navy,dashdotted'


def xfm_current(stem, fault):
    rows = read(stem.relative_to(P).as_posix() + '_current.csv')
    points = [(float(r['hour']) - fault, 100 * float(r['max_current_ratio']))
              for r in rows if r['element'] == 'Transformer.xfm1']
    assert [hour for hour, _ in points] == [i / 4 - fault for i in range(288)]
    return points


def replay(event, case, settings):
    model = Model(settings, 'voltage_only', linear_data(settings))
    stem = DIRECTORY / f'{event}_{case}_voltage_only'
    solution = np.load(stem.with_suffix('.npz'))['x']
    metrics, samples, nodes, _, _, dispatch = evaluate(settings, model, solution)
    archive = read(stem.relative_to(P).as_posix() + '_replay.csv')
    for field in ('max_current_ratio', 'min_energized_voltage', 'max_energized_voltage'):
        assert np.allclose([row[field] for row in samples], [float(row[field]) for row in archive], atol=1e-6)
    assert np.allclose(dispatch['load'].sum(axis=0), [row['care_kw'] for row in samples[::4]])
    fuel = .30 * dispatch['generator'] + 1.5 * dispatch['generator_on']
    assert np.isclose(fuel.sum(), metrics['fuel_litres'])
    energized = {(row['hour'], row['node']): row['voltage_pu'] for row in nodes if row['energized']}
    return dict(load=dispatch['load'], fuel=fuel, generator=dispatch['generator'],
                battery=dispatch['discharge'], voltage=energized, samples=samples,
                current=xfm_current(stem, settings['fault_hour']))


def event_data(event):
    surge, _ = settings_for(event, 'long_reserve')
    base, _ = settings_for(event, 'long_no_surge')
    for key in surge.keys() - {'A'}:
        assert np.array_equal(surge[key], base[key]), key
    fault = surge['fault_hour']
    assert (fault, surge['tie_delay_h'], surge['repair_h']) == (12, 18, 27)
    added = surge['A'] - base['A']
    assert np.all(added >= 0) and not np.any(added[:, :, :fault]) and not np.any(added[:, :, fault+12:])
    cases = {'surge': replay(event, 'long_reserve', surge),
             'base': replay(event, 'long_no_surge', base)}
    for case in cases.values():
        assert np.allclose(case['fuel'][:2, fault:fault+18].sum(axis=1), [320, 320])
    cases['care_current'] = xfm_current(DIRECTORY / f'{event}_long_reserve_coordinated', fault)
    v_surge, v_base = (cases[case]['voltage'] for case in ('surge', 'base'))
    assert v_surge.keys() == v_base.keys()
    cases['max_voltage_change'] = max(abs(voltage - v_base[key]) for key, voltage in v_surge.items())
    cases['arrivals'] = base['A'].sum(axis=(0, 1))
    cases['added'] = added.sum(axis=(0, 1))
    return cases


def coordinates(points, precision=3):
    return ' '.join(f'({x:g},{y:.{precision}f})' for x, y in points)


def step(values):
    return [(hour - 12, value) for hour, value in enumerate([*values, values[-1]])]


def line(style, points, stepwise=True):
    interpolation = ',const plot' if stepwise else ''
    return (rf'\addplot[{style}{interpolation},line width=.76pt] coordinates {{' +
            coordinates(points) + '};')


def annotation(x, y, value, options=''):
    return rf'\node[font=\scriptsize,{options}] at (axis cs:{x},{y}) {{{value}}};'


def resource_table(data):
    rows = [r'\begin{tabular}{llrrr}\toprule',
            r'Storm and case & Site & Diesel (L) & Generation (kWh) & Battery discharge (kWh)\\\midrule']
    for event, title in EVENTS:
        for key, label in (('base', 'No surge'), ('surge', 'Surge')):
            values = data[event][key]
            for site in range(3):
                name = f'{title}, {label}' if site == 0 else ''
                totals = [values[field][site].sum() for field in ('fuel', 'generator', 'battery')]
                rows.append(f'{name} & H{site+1} & ' + ' & '.join(f'{value:.1f}' for value in totals) + r'\\')
    return '\n'.join([*rows, r'\bottomrule\end{tabular}'])


def compact_load_figure(data):
    body = [r'''\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 4,horizontal sep=.84cm,vertical sep=.27cm,
 xlabels at=edge bottom,xticklabels at=edge bottom},
 width=.415\textwidth,height=1.01cm,scale only axis,xmin=-12,xmax=60,
 xtick={-12,0,18,27,42,60},xlabel={Hours from fault},
 tick label style={font=\scriptsize},label style={font=\scriptsize},
 title style={font=\footnotesize\bfseries},ymajorgrids,grid style={gray!20},
 axis line style={black!65}]
''']
    specs = [('Arrivals/h', 0, 40, '0,15,30'), (r'$\Delta$ kW', -40, 85, '-30,0,30,60'),
             (r'XFM1 (\%)', 73, 124, '80,100,120'), ('V (p.u.)', .94, 1.066, '.95,1,1.05')]
    site_styles = ('purple,solid', 'teal,dashed', 'orange!75!black,solid')
    for panel, (label, lo, hi, ticks) in enumerate(specs):
        for event, title in EVENTS:
            d = data[event]
            options = [f'ymin={lo},ymax={hi},ytick={{{ticks}}}']
            if panel == 0:
                options.append(f'title={{{title}}}')
            if event == EVENTS[0][0]:
                options.append(f'ylabel={{{label}}}')
            body.append(r'\nextgroupplot[' + ','.join(options) + ']')
            body.append(rf'\fill[navy!6] (axis cs:0,{lo}) rectangle (axis cs:18,{hi});')
            for hour in (0, 18, 27):
                body.append(rf'\draw[black!55,densely dashed] (axis cs:{hour},{lo})--(axis cs:{hour},{hi});')
            if panel == 0:
                for hour, (base, added) in enumerate(zip(d['arrivals'], d['added']), start=-12):
                    body.append(rf'\fill[teal!43] (axis cs:{hour+.08:g},0) rectangle (axis cs:{hour+.92:g},{base});')
                    if added:
                        body.append(rf'\fill[orange!80!black] (axis cs:{hour+.08:g},{base}) rectangle (axis cs:{hour+.92:g},{base+added});')
                body.append(annotation(58, 34, 'Baseline + added arrivals', 'anchor=east'))
            elif panel == 1:
                delta = d['surge']['load'] - d['base']['load']
                assert lo <= delta.min() <= delta.max() <= hi
                for site, style in enumerate(site_styles):
                    body.append(line(style, step(delta[site])))
                    color = style.split(',')[0]
                    body.append(annotation(34+9*site, 65, f'H{site+1}', f'text={color},anchor=west'))
            elif panel == 2:
                body.append(r'\draw[black!65,densely dotted] (axis cs:-12,100)--(axis cs:60,100);')
                for series, style in ((d['care_current'], CARE), (d['surge']['current'], SURGE),
                                      (d['base']['current'], BASE)):
                    body.append(line(style, [*series, (60, series[-1][1])]))
                peaks = [max(value for _, value in d[case]['current']) for case in ('base', 'surge')]
                hours = [sum(value > 100 + 1e-5 for _, value in d[case]['current']) / 4
                         for case in ('base', 'surge')]
                body.append(annotation(59, 122, rf'{peaks[0]:.1f}\% $\to$ {peaks[1]:.1f}\%\quad'
                                               rf'{hours[0]:g} $\to$ {hours[1]:g} h over',
                                       'anchor=north east,fill=white,inner sep=1pt'))
            else:
                for hospital, style in ((False, 'gray!58,mark size=.29pt'), (True, 'navy,mark size=.39pt')):
                    points = [(hour - 12, value) for (hour, node), value in d['surge']['voltage'].items()
                              if hour.is_integer() and (node.split('.')[0] in ('h1', 'h2', 'h3')) == hospital]
                    body.append(rf'\addplot[only marks,mark=*,{style}] coordinates {{' +
                                coordinates(points, precision=4) + '};')
                for field in ('min_energized_voltage', 'max_energized_voltage'):
                    points = [(row['hour']-12, row[field]) for row in d['base']['samples']]
                    body.append(line(BASE, [*points, (60, points[-1][1])]))
                for limit in (.95, 1.05):
                    body.append(rf'\draw[black!78,densely dotted] (axis cs:-12,{limit})--(axis cs:60,{limit});')
    return '\n'.join([*body, r'\end{groupplot}', r'\end{tikzpicture}'])


def build_load_figure(data=None):
    if data is None:
        data = {event: event_data(event) for event, _ in EVENTS}
    (P / 'No_Surge_Resource_Table.tex').write_text(resource_table(data), encoding='utf-8')
    return '\n'.join([r'\ifdefined\CompactLoadPanels', compact_load_figure(data), r'\else',
                      r'\ifdefined\OmitDieselPanel', detailed_load_figure(data, include_diesel=False),
                      r'\else', detailed_load_figure(data, include_diesel=True), r'\fi\fi'])


def detailed_load_figure(data, *, include_diesel):
    load_max = max(d[case]['load'].max() for d in data.values() for case in ('surge', 'base'))
    load_limit = 60 * ceil((load_max + 8) / 60)
    ticks = ','.join(str(value) for value in range(0, load_limit+1, 60))
    panels = [('Arrivals/h', 0, 42, '0,15,30')]
    panels += [(f'H{i} (kW)', 0, load_limit, ticks) for i in (1, 2, 3)]
    panels += [('Diesel (L)', 0, 880, '0,400,800'), ('XFM1 (\\%)', 70, 128, '80,100,120'),
               ('Voltage (p.u.)', .94, 1.066, '.95,1,1.05')]
    body = [r'''\providecommand{\LoadPanelHeight}{1.20cm}
\providecommand{\LoadPanelGap}{.35cm}
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 7,horizontal sep=.84cm,vertical sep=\LoadPanelGap,
 xlabels at=edge bottom,xticklabels at=edge bottom},
 width=.415\textwidth,height=\LoadPanelHeight,scale only axis,xmin=-12,xmax=60,
 xtick={-12,0,12,18,27,42,60},xlabel={Hours from fault},
 tick label style={font=\scriptsize},label style={font=\scriptsize},
 title style={font=\footnotesize\bfseries},ymajorgrids,grid style={gray!20},
 axis line style={black!65}]
''']
    if not include_diesel:
        body[0] = body[0].replace('2 by 7', '2 by 6')
    for index, (label, lo, hi, ticks) in enumerate(panels):
        if index == 4 and not include_diesel:
            continue
        for event, title in EVENTS:
            d = data[event]
            options = [f'ymin={lo},ymax={hi},ytick={{{ticks}}}']
            if index == 0:
                options.append(f'title={{{title}}}')
            if event == EVENTS[0][0]:
                options.append(f'ylabel={{{label}}}')
            body.append(r'\nextgroupplot[' + ','.join(options) + ']')
            body.append(rf'\fill[navy!6] (axis cs:0,{lo}) rectangle (axis cs:18,{hi});')
            if index == 0:
                for hour, (base, added) in enumerate(zip(d['arrivals'], d['added']), start=-12):
                    body.append(rf'\fill[teal!43] (axis cs:{hour+.08:g},0) rectangle (axis cs:{hour+.92:g},{base});')
                    if added:
                        body.append(rf'\fill[orange!80!black] (axis cs:{hour+.08:g},{base}) rectangle (axis cs:{hour+.92:g},{base+added});')
            if index == 5:
                for hour, current in d['surge']['current']:
                    if current > 100:
                        body.append(rf'\fill[orange!50] (axis cs:{hour:g},100) rectangle (axis cs:{hour+.25:g},{current:.3f});')
            for hour in (0, 18, 27):
                body.append(rf'\draw[black!55,densely dashed] (axis cs:{hour},{lo})--(axis cs:{hour},{hi});')
            if index == 0:
                body.append(annotation(58, 36, 'Baseline + added arrivals', 'anchor=east'))
            elif 1 <= index <= 3:
                site = index - 1
                base, surge = (d[case]['load'][site] for case in ('base', 'surge'))
                body.extend([line(SURGE, step(surge)), line(BASE, step(base))])
                delta = surge - base
                body.append(annotation(59, hi-7, rf'$\Delta$ {delta.min():+.0f} to {delta.max():+.0f} kW',
                                       'anchor=north east,fill=white,inner sep=1pt'))
                if index == 1:
                    body.append(annotation(1, hi-23, 'H1/H2 islanded', 'anchor=west'))
            elif index == 4:
                for case, style in (('surge', SURGE), ('base', BASE)):
                    cumulative = np.r_[0., np.cumsum(d[case]['fuel'].sum(axis=0))]
                    body.append(line(style, [(h-12, value) for h, value in enumerate(cumulative)], stepwise=False))
                totals = [d[case]['fuel'].sum() for case in ('base', 'surge')]
                body.append(annotation(59, 80, rf'{totals[0]:.0f} $\to$ {totals[1]:.0f} L',
                                       'anchor=east,fill=white,inner sep=1pt'))
            elif index == 5:
                body.append(r'\draw[black!65,densely dotted] (axis cs:-12,100)--(axis cs:60,100);')
                for series, style in ((d['care_current'], CARE), (d['surge']['current'], SURGE), (d['base']['current'], BASE)):
                    body.append(line(style, [*series, (60, series[-1][1])]))
                peak_hour, peak = max(d['surge']['current'], key=lambda point: point[1])
                base_peak = max(value for _, value in d['base']['current'])
                body.append(rf'\fill[orange!65!black] (axis cs:{peak_hour:g},{peak:.3f}) circle[radius=1.1pt];')
                body.append(annotation(59, 126, rf'Peak {base_peak:.1f}\% $\to$ {peak:.1f}\%',
                                       'anchor=north east,fill=white,inner sep=1pt'))
            else:
                points = defaultdict(list)
                for (hour, node), voltage in d['surge']['voltage'].items():
                    if hour.is_integer():
                        key = 'hospital' if node.split('.')[0] in ('h1', 'h2', 'h3') else 'feeder'
                        points[key].append((hour - 12, voltage))
                for key, style in [('feeder', 'gray!58,mark size=.32pt'), ('hospital', 'navy,mark size=.44pt')]:
                    body.append(rf'\addplot[only marks,mark=*,{style}] coordinates {{' +
                                coordinates(points[key], precision=4) + '};')
                for field in ('min_energized_voltage', 'max_energized_voltage'):
                    values = [(row['hour']-12, row[field]) for row in d['base']['samples']]
                    body.append(line(BASE, [*values, (60, values[-1][1])]))
                for limit in (.95, 1.05):
                    body.append(rf'\draw[black!78,densely dotted] (axis cs:-12,{limit})--(axis cs:60,{limit});')
                body.append(annotation(59, .999, rf'max $|\Delta V|$ {d["max_voltage_change"]:.3f} p.u.',
                                       'anchor=east,fill=white,inner sep=1pt'))
    body.append(r'\end{groupplot}')
    legend = [(BASE, 'Voltage only, no surge'), (SURGE, 'Voltage only, surge'),
              (CARE, 'Care first, surge (XFM1)')]
    keys = [rf'\tikz[baseline=-.45ex]\draw[{style},line width=.8pt] (0,0)--(.34,0);\,{label}'
            for style, label in legend]
    body.append(r'\node[anchor=south,font=\scriptsize] at ([xshift=.43\textwidth,yshift=.46cm]group c1r1.north west) {' +
                r'\quad '.join(keys) + '};')
    body.append(r'\end{tikzpicture}')
    return '\n'.join(body)
