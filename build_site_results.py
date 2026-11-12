"""One TikZ figure joining hospital outcomes to load, current, and voltage."""

from collections import defaultdict

from build_fault_figures import IST, color, rect
from run_fault_study import Model, np, read
from severe_stress import DIRECTORY, settings_for


EVENTS = (('dec2023', 'December 2023'), ('apr2024', 'April 2024'))
DEMAND = [(-80, '276983'), (0, 'FFFFFF'), (80, 'AD581D')]
CURRENT_ROWS = (('voltage_only', 'Line.h1supply', 'H1 grid'),
                ('voltage_only', 'Line.h2supply', 'H2 grid'),
                ('voltage_only', 'Line.h3supply', 'H3 grid'),
                ('voltage_only', 'Transformer.xfm1', 'XFM1 voltage-only'),
                ('coordinated', 'Transformer.xfm1', 'XFM1 care-first'))
WIDTH = 6.94
DX = WIDTH / 72
ROW = .30


def node(x, y, value, options=''):
    return rf'\node[{options}] at ({x:.4f},{y:.4f}) {{{value}}};'


def row_box(offset, y):
    return rf'\draw[black!40,line width=.25pt] ({offset:g},{y:g}) rectangle ({offset+WIDTH:.4f},{y+ROW:.4f});'


def stage_lines(offset, low, high):
    return [rf'\draw[black!65,densely dashed,line width=.36pt] ({offset+hour*DX:.4f},{low:g})--({offset+hour*DX:.4f},{high:g});'
            for hour in (12, 30, 39)]


def scale(x, y, width, stops, ticks, title, formatter):
    lo, hi = stops[0][0], stops[-1][0]
    parts = 64
    body = [node(x+width/2, y+.22, title, 'anchor=south')]
    for k in range(parts):
        level = lo + (hi-lo)*(k+.5)/parts
        body.append(rect(x+k*width/parts, y, width/parts+.0001, .16, color(level, stops)))
    body.append(rf'\draw[black!50,line width=.3pt] ({x:g},{y:g}) rectangle ({x+width:g},{y+.16:g});')
    for value in ticks:
        tx = x + width*(value-lo)/(hi-lo)
        body.append(rf'\draw[black!75,line width=.35pt] ({tx:.3f},{y:g})--({tx:.3f},{y-.07:g});')
        body.append(node(tx, y-.075, formatter(value), 'anchor=north'))
    return body


def current_rows(event, results):
    data = {}
    for policy in ('voltage_only', 'coordinated'):
        prefix = DIRECTORY.relative_to(DIRECTORY.parent).as_posix()
        records = read(f'{prefix}/{event}_long_reserve_{policy}_current.csv')
        cells = defaultdict(list)
        xf = set()
        for row in records:
            element = row['element']
            hour = float(row['hour'])
            value = float(row['max_current_ratio'])
            cells[element, int(hour)].append(value)
            if element == 'Transformer.xfm1' and value > 1 + 1e-7:
                xf.add(hour)
        assert len(xf)/4 == float(results[event, 'long_reserve', policy]['overload_h'])
        data[policy] = cells
    return data


def outcomes(event, results):
    settings, _ = settings_for(event, 'long_reserve')
    arrivals = settings['A'][:, 0, :].sum(axis=1)
    timely = []
    for policy in ('energy_first', 'coordinated'):
        prefix = DIRECTORY.relative_to(DIRECTORY.parent).as_posix()
        records = read(f'{prefix}/{event}_long_reserve_{policy}_starts.csv')
        counts = [0, 0, 0]
        for row in records:
            if row['care_class'] == '0' and int(row['start_hour']) - int(row['arrival_hour']) <= 2:
                counts[int(row['hospital'][1:])-1] += int(row['patients'])
        metrics = results[event, 'long_reserve', policy]
        assert sum(counts) == int(metrics['high_timely'])
        assert sum(arrivals) == int(metrics['high_arrivals'])
        timely.append(counts)
    lost = []
    for case in ('long100', 'long_reserve'):
        design, _ = settings_for(event, case)
        assert np.array_equal(settings['A'], design['A'])
        model = Model(design, 'coordinated', [])
        archive = np.load(DIRECTORY / f'{event}_{case}_coordinated.npz')['x']
        assert len(archive) == model.n
        site_hours = (1 - np.rint(archive[model.names['plant']])).sum(axis=1).astype(int)
        assert site_hours.sum()*60 == float(results[event, case, 'coordinated']['critical_unserved_kwh'])
        lost.append(site_hours)
    return arrivals, timely, lost


def scatter(body, x, data):
    lo, hi = 5.45, 6.55
    def y(v):
        return lo + (v - .945) * (hi-lo) / .11
    body.append(rf'\fill[navy!5] ({x+12*DX:.4f},{lo}) rectangle ({x+30*DX:.4f},{hi});')
    for value in (.95, 1.05):
        body.append(rf'\draw[black!75,densely dotted] ({x:g},{y(value):.4f})--({x+WIDTH:.4f},{y(value):.4f});')
    body.extend(stage_lines(x, lo, hi))
    points = {True: [], False: []}
    for (hour, phase), value in data['surge']['voltage'].items():
        is_site = phase.split('.')[0] in ('h1', 'h2', 'h3')
        if hour % (1 if is_site else 3) == 0:
            points[is_site].append(rf'({x+hour*DX:.4f},{y(value):.4f}) circle[radius=.39pt]')
    for is_site, style in ((False, 'gray!65'), (True, 'navy')):
        marks = points[is_site]
        for start in range(0, len(marks), 32):
            body.append(rf'\path[fill={style}] ' + ' '.join(marks[start:start+32]) + ';')
    for field in ('min_energized_voltage', 'max_energized_voltage'):
        series = [(s['hour'], s[field]) for s in data['base']['samples'][::4]]
        path = '--'.join(f'({x+t*DX:.4f},{y(v):.4f})' for t, v in series)
        body.append(r'\draw[teal,densely dashed,line width=.36pt] '+path+';')
    body.append(rf'\draw[black!40,line width=.25pt] ({x:g},{lo}) rectangle ({x+WIDTH:.4f},{hi});')
    if x == 0:
        for value in (.95, 1.05):
            body.append(node(-.06, y(value), f'{value:.2f}', 'anchor=east'))


def build_site_figure(results, data):
    body = [r'\begin{tikzpicture}[x=.96cm,y=1cm,every node/.style={font=\scriptsize,inner sep=1.4pt}]']
    headings = [(10.12, '(a) Hospital demand change, surge minus no surge'),
                (8.89, '(b) Current by element'),
                (6.83, '(c) Sampled energized node-phase voltage (p.u.)'),
                (3.96, '(d) Urgent care and essential-service losses by hospital')]
    for y, title in headings:
        body.append(node(-1.56, y, rf'\textbf{{{title}}}', 'anchor=west'))
    body.append(node(8.40+WIDTH, 6.83, 'Hours from fault', 'anchor=east'))
    for idx, (event, title) in enumerate(EVENTS):
        x = idx * 8.40
        body.append(node(x+WIDTH/2, 10.55, rf'\textbf{{{title}}}', r'font=\small'))
        d = data[event]
        delta = d['surge']['load'] - d['base']['load']
        assert abs(delta).max() <= 80
        for i, values in enumerate(delta):
            yy = 9.69-.33*i
            body.append(node(x-.08, yy+.15, f'H{i+1}', 'anchor=east'))
            for hour, value in enumerate(values):
                body.append(rect(x+hour*DX, yy, DX+.0001, ROW, color(value, DEMAND)))
            body.append(row_box(x, yy))
            body.append(node(x+WIDTH-.045, yy+.15, f'{max(values):+.0f} kW',
                             'anchor=east,fill=white,inner sep=.75pt'))
        body.extend(stage_lines(x, 9.03, 9.99))
        current = current_rows(event, results)
        for i, (policy, element, label) in enumerate(CURRENT_ROWS):
            yy = 8.40-.33*i
            if idx == 0:
                body.append(node(x-.08, yy+.15, label, 'anchor=east'))
            for hour in range(72):
                values = current[policy][element, hour]
                value = max(values) if values else None
                body.append(rect(x+hour*DX, yy, DX+.0001, ROW,
                                 'DADDE1' if value is None else color(value, IST)))
                if value is not None and value > 1+1e-7:
                    start = x+hour*DX
                    body.append(rf'\draw[orange!63!black,line width=.38pt] ({start:.4f},{yy:.4f}) rectangle ({start+DX:.4f},{yy+ROW:.4f});')
            body.append(row_box(x, yy))
        body.extend(stage_lines(x, 7.08, 8.70))
        overload = [float(results[event, 'long_reserve', p]['overload_h'])
                    for p in ('voltage_only', 'coordinated')]
        peak = max(value for _, value in d['surge']['current'])
        body.append(node(x+WIDTH, 8.89, rf'XFM1: {peak:.1f}\% peak, {overload[0]:g} $\to$ {overload[1]:g} h over',
                         'anchor=east,text=orange!55!black'))
        scatter(body, x, d)
        for hour in (0, 12, 30, 39, 54, 72):
            body.append(node(x+hour*DX, 5.41, str(hour-12), 'anchor=north'))
        arrivals, timely, lost = outcomes(event, results)
        body.append(node(x+1.68, 3.57, 'On-time urgent share'))
        body.append(node(x+4.70, 3.57, 'Starts / total'))
        body.append(node(x+6.48, 3.57, 'Hours lost'))
        for i, count in enumerate(arrivals):
            yy = 3.16-.54*i
            before, after = [int(a[i]) for a in timely]
            p0, p1 = (x + 3.35*n/count for n in (before, after))
            body.append(node(x-.08, yy, f'H{i+1}', 'anchor=east'))
            body.append(rf'\draw[black!18] ({x:g},{yy:g})--({x+3.35:.4f},{yy:g});')
            body.append(rf'\draw[black!48,line width=.55pt] ({p0:.4f},{yy:g})--({p1:.4f},{yy:g});')
            body.append(rf'\fill[navy] ({p1:.4f},{yy:g}) circle[radius=1.45pt];')
            body.append(rf'\draw[teal,line width=.55pt,fill=white] ({p0:.4f},{yy:g}) circle[radius=1.07pt];')
            body.append(node(x+4.70, yy, rf'{before} $\to$ {after} / {int(count)}'))
            body.append(node(x+6.48, yy, rf'{lost[0][i]} $\to$ {lost[1][i]}'))
        for fraction, label in ((0,'0'),(.5,'50'),(1,r'100\%')):
            body.append(node(x+3.35*fraction, 1.82, label, 'anchor=north'))
        care, energy = [results[event, 'long_reserve', p] for p in ('coordinated', 'energy_first')]
        dk = int(energy['final_energy'])-int(care['final_energy'])
        du = int(energy['high_timely'])-int(care['high_timely'])
        body.append(node(x+WIDTH/2, 1.34, rf'Energy-first trade-off: +{dk} kWh, ${du}$ urgent starts'))
    body.extend(scale(.48, 4.57, 5.98, DEMAND, [-80, -40, 0, 40, 80],
                      'Demand change (kW)', lambda v: str(v)))
    body.extend(scale(8.88, 4.57, 5.98, IST, [0, .4, .8, 1, 1.2],
                      r'Current / rating (\%)', lambda v: str(round(100*v))))
    body.append(node(7.67, .94, r'\textcolor{teal}{$\circ$ Energy first}\quad\textcolor{navy}{$\bullet$ Care first}\qquad Essential hours: 160 L $\to$ 320 L'))
    body.append(r'\end{tikzpicture}')
    return '\n'.join(body)
