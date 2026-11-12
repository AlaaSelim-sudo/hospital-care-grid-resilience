"""Build the severe-case comparison table and fuel-sizing TikZ plot."""

import math
from collections import Counter, defaultdict

from build_fault_figures import ELEMENTS, IST, color, colorbar, label, rect
from build_load_figure import annotation, build_load_figure, event_data, line
from build_site_results import build_site_figure
from build_care_diesel_figure import diesel_figure
from compound_stress import EVENTS
from run_fault_study import P, read
from severe_stress import CASES, DIRECTORY, POLICIES, settings_for


URGENT_DEADLINE_H = 2


TRACE_POLICIES = (
    ('coordinated', 'navy,solid', 'Care first'),
    ('energy_first', 'teal,dashed', 'Energy first'),
    ('voltage_only', 'orange,densely dotted', 'Voltage only'),
)

DEADLINE_POLICIES = TRACE_POLICIES[:2]
HOSPITAL_GRID_LINES = {'Line.h1supply': 'H1 grid', 'Line.h2supply': 'H2 grid',
                       'Line.h3supply': 'H3 grid'}


def results():
    rows = {(r['event'], r['case'], r['policy']): r
            for name in ('Severe_Stress_Results.csv', 'No_Surge_Results.csv') for r in read(name)}
    assert rows.keys() == {(e, c, p) for e in EVENTS for c in CASES for p in POLICIES[c]}
    reserves = {(r['event'], int(r['isolation_h'])): r for r in read('Severe_Reserve_Results.csv')}
    assert reserves.keys() == {(e, t) for e in EVENTS for t in (15, 18)}
    return rows, reserves


def table(rows):
    lines = [r'\begin{tabular}{llrrrrrr}\toprule',
             r'Event & Policy & Island (h) & Surge (\%) & Fuel (L) & Deficit & Overcurrent (h) & Urgent\\\midrule']
    labels = {'coordinated': 'Care', 'energy_first': 'Energy', 'voltage_only': 'Voltage'}
    for event, label in zip(EVENTS, ('Dec.', 'Apr.')):
        for case in CASES:
            if case == 'long_no_surge':
                continue
            for policy in POLICIES[case]:
                row = rows[event, case, policy]
                essential = row['critical_unserved_kwh']
                urgent = '---'
                if (row['AC_voltage_thermal_feasible'] == 'True' and row['essential_certified'] == 'True'
                        and (policy != 'energy_first' or row['energy_certified'] == 'True')):
                    upper = row.get('urgent_upper_bound')
                    urgent = row['high_timely'] + ('/' + upper if upper else '/---')
                policy_label = 'Full pickup' if case.startswith('long_full') else labels[policy]
                values = (label, policy_label, row['isolation_h'], row['surge_pct'], row['fuel_per_island_l'],
                          f'{float(essential):.0f}' if essential else '---', row['overload_h'] or '---', urgent)
                lines.append(' & '.join(values) + r'\\')
                label = ''
        lines.append(r'\addlinespace')
    return '\n'.join(lines + [r'\bottomrule\end{tabular}'])


def fuel_plot(reserves):
    points = []
    for t in (15, 18):
        r = reserves['dec2023', t]
        other = reserves['apr2024', t]
        assert r['AC_pass'] == other['AC_pass'] == 'True'
        stock = float(r['capacity_per_island_l'])
        assert math.isclose(stock, float(other['capacity_per_island_l']), abs_tol=1e-5)
        points.append(f'({t},{stock:.3f})')
    physical = []
    for t in (15, 18):
        generation = 60 * t - .95 * (200 - 20)
        litres = .30 * generation + 1.5 * math.ceil(generation / 100)
        physical.append(f'({t},{litres:.3f})')
    return r'''\begin{tikzpicture}
\begin{axis}[width=.74\linewidth,height=5.4cm,xtick={15,18},grid=major,
 xlabel={Hours of isolation},ylabel={Fuel per islanded site (L)},
 legend style={at={(.5,1.04)},anchor=south,legend columns=2,draw=none,font=\small}]
\addplot[navy,thick,mark=*] coordinates {''' + ' '.join(points) + r'''};
\addlegendentry{Guarded model, both storms}
\addplot[black!65,dashed,mark=triangle*] coordinates {''' + ' '.join(physical) + r'''};
\addlegendentry{Physical lower bound}
\end{axis}\end{tikzpicture}'''


def timely_urgent_by_hour(prefix, group):
    counts = Counter()
    for row in read(prefix + '_starts.csv'):
        arrival = int(row['arrival_hour'])
        start = int(row['start_hour'])
        if row['care_class'] == '0' and start - arrival <= URGENT_DEADLINE_H:
            counts[int(row[group])] += int(row['patients'])
    return counts


def load_traces(rows):
    traces = {}
    for event in EVENTS:
        for policy, _, _ in TRACE_POLICIES:
            prefix = DIRECTORY.relative_to(P).as_posix() + f'/{event}_long_reserve_{policy}'
            samples = read(prefix + '_replay.csv')
            metrics = rows[event, 'long_reserve', policy]
            fault = next(float(r['hour']) for r in samples if r['stage'] == 'island')
            tie = next(float(r['hour']) for r in samples if r['stage'] == 'backfeed')
            repair = next(float(r['hour']) for r in samples if float(r['hour']) > tie and r['stage'] == 'normal')
            assert (fault, tie - fault, repair - fault) == (12, 18, 27)
            assert [float(r['hour']) for r in samples] == [i / 4 for i in range(288)]
            x = [float(r['hour']) - fault for r in samples]
            low_v = [float(r['min_energized_voltage']) for r in samples]
            high_v = [float(r['max_energized_voltage']) for r in samples]
            current = [100 * float(r['max_current_ratio']) for r in samples]
            assert all(v >= .95 - 1e-7 for v in low_v) and all(v <= 1.05 + 1e-7 for v in high_v)
            assert float(metrics['voltage_violation_h']) == 0

            over = [t for t, row, ratio in zip(x, samples, current)
                    if int(row['overloaded_current_entries']) > 0 and ratio > 100]
            assert len(over) / 4 == float(metrics['overload_h'])
            intervals = []
            for t in over:
                if intervals and intervals[-1][1] == t:
                    intervals[-1] = (intervals[-1][0], t + .25)
                else:
                    intervals.append((t, t + .25))

            starts = timely_urgent_by_hour(prefix, 'start_hour')
            cumulative = [(x[0], 0)]
            for start, patients in sorted(starts.items()):
                cumulative.append((start - fault, cumulative[-1][1] + patients))
            assert cumulative[-1][1] == int(metrics['high_timely'])
            end = x[-1] + .25
            traces[event, policy] = dict(voltage=list(zip(x + [end], low_v + [low_v[-1]])),
                                         current=list(zip(x + [end], current + [current[-1]])),
                                         urgent=cumulative + [(end, cumulative[-1][1])], over=intervals)
    return traces


def time_plot(rows):
    traces = load_traces(rows)
    body = [r'''\begin{tikzpicture}
\begin{groupplot}[
 group style={group size=2 by 3,horizontal sep=.95cm,vertical sep=.21cm,
              xlabels at=edge bottom,xticklabels at=edge bottom},
 width=.405\textwidth,height=2.23cm,scale only axis,xmin=-12,xmax=60,
 xtick={-12,0,18,27,42,60},xlabel={Hours from fault},
 tick label style={font=\scriptsize},label style={font=\scriptsize},
 title style={font=\footnotesize\bfseries},ymajorgrids,grid style={gray!22},
 axis line style={black!70},unbounded coords=jump,
 after end axis/.code={\foreach \stagepos in {.166666667,.416666667,.541666667}{%
  \draw[gray,densely dashed,thin] (rel axis cs:\stagepos,0)--(rel axis cs:\stagepos,1);}}]
''']
    panels = [
        ('voltage', 'Min. voltage (p.u.)', '.945', '1.005', '{.95,.97,.99}', '.95'),
        ('current', 'Max. current (\\%)', '78', '115', '{80,90,100,110}', '100'),
        ('urgent', 'Timely urgent starts', '0', '112', '{0,50,100}', None),
    ]
    for field, axis_label, low, high, ticks, limit in panels:
        for event, title in zip(EVENTS, ('December 2023', 'April 2024')):
            options = [f'ymin={low}', f'ymax={high}', f'ytick={ticks}']
            if event == EVENTS[0]:
                options.append(f'ylabel={{{axis_label}}}')
            if field == 'voltage':
                options.append(f'title={{{title}}}')
            body.append(r'\nextgroupplot[' + ','.join(options) + ']')
            if field == 'current':
                for start, stop in traces[event, 'voltage_only']['over']:
                    body.append(rf'\fill[orange!19] (axis cs:{start:g},78) rectangle (axis cs:{stop:g},115);')
            if limit is not None:
                body.append(rf'\addplot[black!65,densely dashed,thin] coordinates {{(-12,{limit}) (60,{limit})}};')
            for policy, style, _ in TRACE_POLICIES:
                points = ' '.join(f'({t:g},{value:.5f})' for t, value in traces[event, policy][field])
                body.append(rf'\addplot[{style},line width=.7pt,const plot] coordinates {{{points}}};')
    body.append(r'\end{groupplot}')
    legend = [rf'\tikz[baseline=-.45ex]\draw[{style},line width=.8pt] (0,0)--(.33,0);\,{label}'
              for _, style, label in TRACE_POLICIES]
    body.append(r'''\node[anchor=south,font=\scriptsize] at
 ([xshift=.455\textwidth,yshift=.46cm]group c1r1.north west) {''' + r'\quad '.join(legend) + '};')
    body.append(r'\end{tikzpicture}')
    return '\n'.join(body)


def current_heatmap(rows):
    event, case = 'apr2024', 'long_reserve'
    figure = [r'\begin{tikzpicture}[x=.092cm,y=.228cm,every node/.style={font=\scriptsize}]']
    for offset, policy, title in [(0, 'voltage_only', 'Voltage only'), (94, 'coordinated', 'Care first')]:
        prefix = DIRECTORY.relative_to(P).as_posix() + f'/{event}_{case}_{policy}_current.csv'
        cells = defaultdict(list)
        overloaded = set()
        for sample in read(prefix):
            t = float(sample['hour'])
            element = sample['element']
            assert element in ELEMENTS and 0 <= t < 72
            ratio = float(sample['max_current_ratio'])
            assert IST[0][0] <= ratio <= IST[-1][0], (policy, t, element, ratio)
            cells[element, int(t)].append(ratio)
            if ratio > 1 + 1e-7:
                overloaded.add(t)
                assert element == 'Transformer.xfm1'
        assert len(overloaded) / 4 == float(rows[event, case, policy]['overload_h'])
        figure.append(rf'\node[font=\small\bfseries] at ({offset+36},22.4) {{{title}}};')
        for index, element in enumerate(ELEMENTS):
            y = len(ELEMENTS) - 1 - index
            style = 'anchor=east,font=\\scriptsize\\bfseries' if element == 'Transformer.xfm1' else 'anchor=east'
            row_label = HOSPITAL_GRID_LINES.get(element, label(element))
            figure.append(rf'\node[{style}] at ({offset-1},{y+.5}) {{{row_label}}};')
            if element in ('Line.h1supply', 'Line.h2supply'):
                assert all(not cells[element, h] for h in range(12, 30))
                assert all(max(cells[element, h]) < .25
                           for h in range(72) if not 12 <= h < 30)
            for h in range(72):
                values = cells[element, h]
                peak = max(values) if values else None
                fill = color(peak, IST) if peak is not None else 'DADDE1'
                figure.append(rect(offset + h, y, 1, 1, fill))
                if peak is not None and peak > 1 + 1e-7:
                    figure.append(rf'\draw[black!35!orange,line width=.36pt] ({offset+h}, {y}) rectangle ({offset+h+1}, {y+1});')
        figure.append(rf'\draw[black!50] ({offset},0) rectangle ({offset+72},21);')
        for elapsed in (0, 12, 24, 30, 39, 54, 72):
            figure.append(rf'\node[anchor=north] at ({offset+elapsed},-.2) {{{elapsed-12}}};')
        for elapsed in (12, 30, 39):
            figure.append(rf'\draw[black!65,densely dashed] ({offset+elapsed},0)--({offset+elapsed},21);')
        figure.append(rf'\node[anchor=north] at ({offset+36},-1.65) {{Hours from fault}};')
    figure.append(colorbar(18, -5, 130, IST, [0, .4, .8, 1, 1.2],
                           r'Peak phase current / rating (\%)', lambda x: f'{100*x:.0f}'))
    figure.append(r'\end{tikzpicture}')
    return '\n'.join(figure)


def deadline_curve(start_rows, settings, metrics):
    fault = settings['fault_hour']
    arrivals = settings['A'][:, 0, :].sum(axis=0)
    timely = Counter()
    for row in start_rows:
        arrival = int(row['arrival_hour'])
        if row['care_class'] == '0' and int(row['start_hour']) - arrival <= URGENT_DEADLINE_H:
            timely[arrival] += int(row['patients'])
    assert metrics['AC_voltage_thermal_feasible'] == 'True'
    assert int(arrivals.sum()) == int(metrics['high_arrivals'])
    assert sum(timely.values()) == int(metrics['high_timely'])
    points = [(-fault, 0)]
    for hour, arrivals_in_hour in enumerate(arrivals):
        missed = int(arrivals_in_hour) - timely[hour]
        assert missed >= 0
        if missed:
            points.append((hour + URGENT_DEADLINE_H - fault, points[-1][1] + missed))
    assert points[-1][0] < 42
    return points + [(42, points[-1][1])]



def essential_service_curve(event, case, settings, metrics):
    prefix = DIRECTORY.relative_to(P).as_posix()
    trace = read(f'{prefix}/{event}_{case}_coordinated_replay.csv')
    site_count = settings['A'].shape[0]
    samples = [round(site_count * float(row['energized_service_fraction'])) for row in trace]
    assert all(len(set(samples[i:i+4])) == 1 for i in range(0, len(samples), 4))
    served = samples[::4]
    assert 60 * sum(site_count - value for value in served) == float(metrics['critical_unserved_kwh'])
    return [(hour - settings['fault_hour'], value) for hour, value in enumerate(served)
            if hour - settings['fault_hour'] <= 42]


def focused_resilience_plot(rows, deadlines, essential, *, include_energy=True):
    body = [r'''\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 3,horizontal sep=.88cm,vertical sep=.39cm},
 width=.405\textwidth,height=1.12cm,scale only axis,ymajorgrids,grid style={gray!20},
 tick label style={font=\scriptsize},label style={font=\scriptsize},
 title style={font=\scriptsize\bfseries},axis line style={black!65}]
''']
    if not include_energy:
        body[0] = body[0].replace('2 by 3', '2 by 2')
    for panel in range(3 if include_energy else 2):
        for event, title in zip(EVENTS, ('December 2023', 'April 2024')):
            if panel == 0:
                options = ['xmin=-12,xmax=48,ymin=0,ymax=65,ytick={0,20,40},xtick={-12,0,18,27,42}',
                           'xticklabels=\\empty', f'title={{{title}}}']
            elif panel == 1:
                options = ['xmin=-12,xmax=48,ymin=.7,ymax=3.25,ytick={1,2,3},xtick={-12,0,18,27,42}']
            else:
                options = ['xmin=4160,xmax=4390,ymin=80,ymax=108,ytick={85,100},xtick={4200,4300,4380}',
                           r'xticklabel style={/pgf/number format/1000 sep={}}']
            if event == EVENTS[0]:
                options.append('ylabel={' + ('Misses', 'Sites served', 'On time')[panel] + '}')
            body.append(r'\nextgroupplot[' + ','.join(options) + ']')
            if panel < 2:
                lo, hi = ((0, 65), (.7, 3.25))[panel]
                for hour in (0, 18, 27):
                    body.append(rf'\draw[black!50,densely dashed] (axis cs:{hour},{lo})--(axis cs:{hour},{hi});')
            if panel == 0:
                for policy, style, _ in DEADLINE_POLICIES:
                    points = deadlines[event, policy]
                    body.append(line(style, points))
                    last = points[-1][1]
                    color_name = style.split(',')[0]
                    body.append(annotation(42, last, str(last), f'anchor=west,text={color_name},inner sep=1pt'))
                body.append(annotation(-10, 57, r'\textcolor{navy}{Care first}\quad\textcolor{teal}{Energy first}',
                                       'anchor=west,fill=white,inner sep=1pt'))
            elif panel == 1:
                for hour, served in essential[event, 'long100']:
                    if served < 3:
                        body.append(rf'\fill[orange!24] (axis cs:{hour},{served}) rectangle (axis cs:{hour+1},3);')
                body.append(line('orange!60!black,dashed', essential[event, 'long100']))
                body.append(line('navy,solid', essential[event, 'long_reserve']))
                gap = float(rows[event, 'long100', 'coordinated']['critical_unserved_kwh'])
                gap -= float(rows[event, 'long_reserve', 'coordinated']['critical_unserved_kwh'])
                body.append(annotation(46, 1.48, f'+{gap / 60:g} hospital-hours',
                                       'anchor=east,fill=white,inner sep=1pt'))
                body.append(annotation(46, 2.28, r'\textcolor{orange!60!black}{160 L}\enspace\textcolor{navy}{320 L}',
                                       'anchor=east,fill=white,inner sep=1pt'))
            else:
                care, energy = [rows[event, 'long_reserve', policy]
                                for policy in ('coordinated', 'energy_first')]
                points = [(int(r['final_energy']), int(r['high_timely'])) for r in (care, energy)]
                (x0, y0), (x1, y1) = points
                body.append(line('black!45', points, stepwise=False))
                for (x, y), color_name in zip(points, ('navy', 'teal')):
                    body.append(rf'\fill[{color_name}] (axis cs:{x},{y}) circle[radius=1.25pt];')
                body.append(annotation(4170, 86, rf'+{x1-x0} kWh, ${y1-y0}$ starts',
                                       'anchor=west,fill=white,inner sep=1pt'))
    return '\n'.join([*body, r'\end{groupplot}', r'\end{tikzpicture}'])


def resilience_plot(rows):
    lines = [r'''\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 3,horizontal sep=.88cm,vertical sep=.43cm,
 xlabels at=edge bottom,xticklabels at=edge bottom},
 width=.405\textwidth,height=1.22cm,scale only axis,
 xmin=-12,xmax=47,xtick={-12,0,12,18,27,42},
 xlabel={Hours from fault},ymajorgrids,grid style={gray!20},
 tick label style={font=\scriptsize},label style={font=\scriptsize},
 title style={font=\scriptsize\bfseries},axis line style={black!65},
 legend style={font=\scriptsize,draw=none,fill=white,inner sep=1pt,
               legend columns=2,/tikz/every even column/.append style={column sep=2pt}}]
''']
    arrivals, deadlines, essential = {}, {}, {}
    for event in EVENTS:
        settings, _ = settings_for(event, 'long_reserve')
        low, _ = settings_for(event, 'long100')
        assert (settings['A'] == low['A']).all() and settings['stages'] == low['stages']
        assert settings['backfeed_background_fraction'] == low['backfeed_background_fraction']
        care = rows[event, 'long_reserve', 'coordinated']
        energy = rows[event, 'long_reserve', 'energy_first']
        for field in ('arrivals', 'high_arrivals', 'fuel_per_island_l', 'isolation_h',
                      'repair_h', 'backfeed_native_fraction', 'final_essential'):
            assert care[field] == energy[field]
        assert care['fuel_per_island_l'] == '320'
        urgent, other = settings['A'].sum(axis=0)
        assert int(urgent.sum()) == int(care['high_arrivals'])
        assert int(urgent.sum() + other.sum()) == int(care['arrivals'])
        arrivals[event] = [(hour - settings['fault_hour'], int(u), int(o))
                           for hour, (u, o) in enumerate(zip(urgent, other))]
        for policy, _, _ in DEADLINE_POLICIES:
            prefix = DIRECTORY.relative_to(P).as_posix()
            starts = read(f'{prefix}/{event}_long_reserve_{policy}_starts.csv')
            deadlines[event, policy] = deadline_curve(starts, settings, rows[event, 'long_reserve', policy])
        for case in ('long_reserve', 'long100'):
            essential[event, case] = essential_service_curve(
                event, case, settings, rows[event, case, 'coordinated'])

    specifications = [('Arrivals/h', 40, '0,15,30'), ('Urgent misses', 65, '0,20,40'),
                      ('Sites served', 3.25, '1,2,3')]
    fuel_styles = [('long100', 'orange!60!black,dashed', '160 L'), ('long_reserve', 'navy,solid', '320 L')]
    for index, (label, ymax, ticks) in enumerate(specifications):
        for event, title in zip(EVENTS, ('December 2023', 'April 2024')):
            options = [f'ymin=0,ymax={ymax},ytick={{{ticks}}}']
            if index == 0:
                options.append(f'title={{{title}}}')
            if event == EVENTS[0]:
                options.append(f'ylabel={{{label}}}')
            legend_x, legend_y, anchor = [('.98', '.97', 'north east'),
                                          ('.02', '.97', 'north west'),
                                          ('.98', '.03', 'south east')][index]
            options.append(rf'legend style={{at={{({legend_x},{legend_y})}},anchor={anchor}}}')
            lines.append(r'\nextgroupplot[' + ','.join(options) + ']')

            if index == 0:
                lines.append(r'\fill[gray!14] (axis cs:0,0) rectangle (axis cs:12,40);')
                for hour, urgent, other in arrivals[event]:
                    assert urgent + other < ymax
                    if urgent:
                        lines.append(rf'\fill[navy!82] (axis cs:{hour+.08:g},0) rectangle (axis cs:{hour+.92:g},{urgent});')
                    if other:
                        lines.append(rf'\fill[teal!37] (axis cs:{hour+.08:g},{urgent}) rectangle (axis cs:{hour+.92:g},{urgent+other});')
                if event == EVENTS[1]:
                    lines.extend([r'\addlegendimage{area legend,fill=navy!82,draw=none}',
                                  r'\addlegendentry{Urgent}',
                                  r'\addlegendimage{area legend,fill=teal!37,draw=none}',
                                  r'\addlegendentry{Other}'])

            if index == 2:
                for hour, served in essential[event, 'long100']:
                    if served < 3:
                        lines.append(rf'\fill[orange!20] (axis cs:{hour},{served}) rectangle (axis cs:{hour+1},3);')
            for hour in (0, 18, 27):
                lines.append(rf'\draw[black!50,densely dashed] (axis cs:{hour},0)--(axis cs:{hour},{ymax});')
            if index == 0:
                continue
            for policy, style, legend in DEADLINE_POLICIES if index == 1 else fuel_styles:
                points = deadlines[event, policy] if index == 1 else essential[event, policy]
                coordinates = ' '.join(f'({hour},{value:g})' for hour, value in points)
                lines.append(rf'\addplot[{style},line width=.85pt,const plot] coordinates {{{coordinates}}};')
                if index == 1:
                    color_name = style.split(',')[0]
                    last = points[-1][1]
                    lines.append(rf'\node[anchor=west,inner sep=1pt,font=\scriptsize\bfseries,text={color_name}] at (axis cs:42,{last}) {{{last}}};')
                if event == EVENTS[1]:
                    lines.append(rf'\addlegendentry{{{legend}}}')
    lines.extend([r'\end{groupplot}', r'\end{tikzpicture}'])
    return '\n'.join([r'\ifdefined\ConferenceResiliencePanels',
                      focused_resilience_plot(rows, deadlines, essential, include_energy=False), r'\else',
                      r'\ifdefined\FocusResiliencePanels', focused_resilience_plot(rows, deadlines, essential),
                      r'\else', *lines, r'\fi\fi'])

def build():
    rows, reserves = results()
    (P / 'Severe_Stress_Table.tex').write_text(table(rows), encoding='utf-8')
    (P / 'figures' / 'severe_fuel_requirements.tikz').write_text(fuel_plot(reserves), encoding='utf-8')
    (P / 'figures' / 'severe_time_series.tikz').write_text(time_plot(rows), encoding='utf-8')
    (P / 'figures' / 'severe_current_heatmap.tikz').write_text(current_heatmap(rows), encoding='utf-8')
    (P / 'figures' / 'severe_urgent_delays.tikz').write_text(resilience_plot(rows), encoding='utf-8')
    load_data = {event: event_data(event) for event in EVENTS}
    (P / 'figures' / 'severe_load_stress.tikz').write_text(build_load_figure(load_data), encoding='utf-8')
    (P / 'figures' / 'site_resilience_composite.tikz').write_text(build_site_figure(rows, load_data), encoding='utf-8')
    diesel_figure(deadline_curve)
    return rows, reserves


def manuscript_values(rows, reserves):
    values = {}
    for event in EVENTS:
        for case in ('mid50', 'mid100', 'long100', 'long_reserve'):
            r = rows[event, case, 'coordinated']
            assert r['AC_voltage_thermal_feasible'] == r['essential_certified'] == 'True'
            assert abs(float(r['critical_unserved_kwh']) - float(r['physical_deficit_floor_kwh'])) < 1e-5
        for hours in (15, 18):
            r = reserves[event, hours]
            assert r['AC_pass'] == 'True' and float(r['solver_gap']) < 1e-5
        assert rows[event, 'long_full', 'coordinated']['essential_status'] == '8'
        assert rows[event, 'long_full_reserve', 'coordinated']['essential_status'] == '8'
        for policy, token in (('coordinated', 'CARE'), ('energy_first', 'ENERGY')):
            row = rows[event, 'long_reserve', policy]
            values[f'MISSED_{event.upper()}_{token}'] = str(
                int(row['high_arrivals']) - int(row['high_timely']))

    values['DEADLINE_GAIN'] = ' and '.join(
        str(int(values[f'MISSED_{event.upper()}_ENERGY']) - int(values[f'MISSED_{event.upper()}_CARE']))
        for event in EVENTS)

    values['SEVERE_FLOOR'] = ' and '.join(rows['dec2023', c, 'coordinated']['physical_deficit_floor_kwh']
                                    for c in ('mid100', 'long100'))
    stocks = [[f"{float(reserves[e, t]['capacity_per_island_l']):.1f}" for t in (15, 18)] for e in EVENTS]
    assert stocks[0] == stocks[1]
    values['SEVERE_STOCK'] = ' and '.join(stocks[0])
    lower = [rows[e, 'mid50', 'coordinated'] for e in EVENTS]
    upper = [rows[e, 'mid100', 'coordinated'] for e in EVENTS]
    values['MID_LOWER'] = ' and '.join(f"{100*int(r['high_timely'])/int(r['high_arrivals']):.1f}\\%" for r in lower)
    values['MID_UPPER'] = ' and '.join(f"{100*int(r['urgent_upper_bound'])/int(r['high_arrivals']):.1f}\\%" for r in upper)
    values['SEVERE_URGENT_GAIN'] = ' and '.join(
        str(int(rows[e, 'long_reserve', 'coordinated']['high_timely']) -
            int(rows[e, 'long100', 'coordinated']['urgent_upper_bound'])) for e in EVENTS)
    values['SEVERE_OVERLOAD'] = ' and '.join(
        f"{float(rows[e, 'long_reserve', 'voltage_only']['overload_h']):g}" for e in EVENTS)
    values['NO_SURGE_OVERLOAD'] = ' and '.join(
        f"{float(rows[e, 'long_no_surge', 'voltage_only']['overload_h']):g}" for e in EVENTS)
    tradeoff = {r['event']: r for r in read('Severe_Tradeoff_Results.csv')}
    for event in EVENTS:
        r = tradeoff[event]
        assert r['proven'] == 'True'
        assert int(r['upper_at_urgent_max_kwh']) < int(r['energy_candidate_kwh'])
    for token, field in [('SEVERE_URGENT_MAX', 'urgent_maximum'), ('SEVERE_ENERGY_CAP', 'upper_at_urgent_max_kwh'),
                         ('SEVERE_ENERGY_WITNESS', 'energy_candidate_kwh')]:
        values[token] = ' and '.join(f'{int(tradeoff[e][field]):,}' for e in EVENTS)
    for event in EVENTS:
        value = int(rows[event, 'long_reserve', 'coordinated']['final_energy'])
        values[f'CARE_ENERGY_{event.upper()}'] = f'{value:,}'
    energy_runs = [rows[event, 'long_reserve', 'energy_first'] for event in EVENTS]
    def interval(run, stage):
        value, upper = (int(run[f'{stage}_{key}']) for key in ('value', 'upper_bound'))
        assert value <= upper
        return str(value) if value == upper else f'{value}--{upper}'
    energy_ranges = ' and '.join(interval(run, 'energy') for run in energy_runs)
    urgent_ranges = ' and '.join(interval(run, 'urgent') for run in energy_runs)
    if all(run['energy_certified'] == run['urgent_certified'] == 'True' for run in energy_runs):
        values['ENERGY_PROOF_NOTE'] = 'Both energy optima and urgent-care tie-breaks were certified.'
    elif all(run['energy_certified'] == 'True' for run in energy_runs):
        values['ENERGY_PROOF_NOTE'] = (
            'Both energy optima were certified. Their urgent tie-breaks give '
            f'{urgent_ranges} on-time starts, respectively.')
    else:
        values['ENERGY_PROOF_NOTE'] = (
            f'The energy optima lie within {energy_ranges} kWh. At the plotted energy levels, '
            f'the on-time urgent optima lie within {urgent_ranges}, in storm order.')
    return values
