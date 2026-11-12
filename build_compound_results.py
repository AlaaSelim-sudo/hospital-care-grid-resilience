"""Generate the compound-case table and native TikZ traces from archived schedules."""

from compound_stress import DIRECTORY, EVENTS, POLICIES, design
from run_fault_study import P, np, read


def results():
    rows = {(r['event'], r['pattern'], r['policy']): r for r in read('Compound_Stress_Results.csv')}
    expected = {(e, s, p) for e in EVENTS for s in ('synchronized', 'extended') for p in POLICIES}
    expected |= {(e, s, 'coordinated') for e in EVENTS for s in ('staggered', 'reserve')}
    assert rows.keys() == expected
    return rows


def april_tradeoff(rows):
    energy = rows['apr2024', 'extended', 'energy_first']
    care = rows['apr2024', 'extended', 'coordinated']
    assert energy['essential_certified'] == care['essential_certified'] == care['urgent_certified'] == 'True'
    assert energy['final_essential'] == care['final_essential']
    assert energy['AC_voltage_thermal_feasible'] == care['AC_voltage_thermal_feasible'] == 'True'
    urgent = int(care['urgent_value'])
    assert int(care['final_urgent']) == urgent
    capped_energy = int(care['energy_upper_bound'])
    feasible_energy = int(energy['final_energy'])
    assert capped_energy < feasible_energy <= int(energy['energy_upper_bound'])
    return urgent, capped_energy, feasible_energy


def extended_table(rows):
    table = [r'\begin{tabular}{llrrrr}\toprule',
             r'Event & Policy & Urgent & Deficit & $V$ (h) & $I$ (h)\\\midrule']
    for event, label in zip(EVENTS, ('Dec.', 'Apr.')):
        for pattern, policy, name in [('extended', 'autonomous', 'Local'),
                                      ('extended', 'voltage_only', 'Voltage'),
                                      ('extended', 'energy_first', 'Energy'),
                                      ('extended', 'coordinated', 'Care'),
                                      ('reserve', 'coordinated', 'Care + fuel')]:
            row = rows[event, pattern, policy]
            urgent = row['high_timely'] if row['AC_voltage_thermal_feasible'] == 'True' else '--'
            if policy == 'energy_first' and row['energy_certified'] != 'True':
                urgent = r'\textit{open}'
            cells = (label, name, urgent, f"{float(row['critical_unserved_kwh']):.0f}",
                     f"{float(row['voltage_violation_h']):.2f}", f"{float(row['overload_h']):.2f}")
            table.append(' & '.join(cells) + r'\\')
            label = ''
    return '\n'.join(table + [r'\bottomrule\end{tabular}'])


def coordinates(points):
    return ' '.join(f'({x:g},{y:.4f})' for x, y in points)


def traces(event):
    plots = [('extended', 'voltage_only', 'orange,densely dotted', 'Voltage only'),
             ('extended', 'coordinated', 'navy', 'Care first'),
             ('reserve', 'coordinated', 'teal,dashdotted', 'Care + fuel')]
    electrical, clinical = [], []
    for pattern, policy, style, label in plots:
        stem = DIRECTORY.relative_to(P).as_posix() + f'/{event}_{pattern}_{policy}'
        replay = read(stem + '_replay.csv')
        current = [(float(r['hour']), 100 * float(r['max_current_ratio'])) for r in replay
                   if 8 <= float(r['hour']) <= 40]
        electrical.append(f'\\addplot[{style},thick,const plot] coordinates {{{coordinates(current)}}};')
        electrical.append(f'\\addlegendentry{{{label}}}')
        if policy == 'voltage_only':
            continue
        if any(int(r['violating_node_phases']) or int(r['overloaded_current_entries']) for r in replay):
            continue
        settings, _, _ = design(event, pattern)
        arrivals = settings['A'][:, 0, :].sum(axis=0)
        timely = np.zeros(arrivals.size, dtype=int)
        for row in read(stem + '_starts.csv'):
            arrival = int(row['arrival_hour'])
            if int(row['care_class']) == 0 and int(row['start_hour']) - arrival <= 2:
                timely[arrival] += int(row['patients'])
        missed = [(i + 2, y) for i, y in enumerate(np.cumsum(arrivals - timely)) if 8 <= i + 2 <= 40]
        clinical.append(f'\\addplot[{style},thick,const plot] coordinates {{{coordinates(missed)}}};')
    return r'''\begin{tikzpicture}
\begin{groupplot}[group style={group size=1 by 2,vertical sep=1.18cm},width=.91\linewidth,height=4.45cm,
 xmin=8,xmax=40,xtick={8,12,16,20,24,28,33,36,40},tick label style={font=\small},
 label style={font=\small},grid=major,clip=false,
 after end axis/.code={\foreach \position in {.125,.5,.78125}{
 \draw[black!40,densely dashed] (rel axis cs:\position,0)--(rel axis cs:\position,1);}}]
\nextgroupplot[ylabel={Max. current (\%)},legend style={at={(.5,1.13)},anchor=south,legend columns=3,font=\footnotesize,draw=none}]
\addplot[black!65,densely dashed,forget plot] coordinates {(8,100) (40,100)};
''' + '\n'.join(electrical) + r'''
\nextgroupplot[ylabel={Cumulative urgent misses},xlabel={Elapsed hour},ymin=0]
''' + '\n'.join(clinical) + r'''
\end{groupplot}\end{tikzpicture}'''


def main():
    rows = results()
    (P / 'Compound_Stress_Table.tex').write_text(extended_table(rows), encoding='utf-8')
    (P / 'figures' / 'compound_stress_april.tikz').write_text(traces('apr2024'), encoding='utf-8')


if __name__ == '__main__':
    main()
