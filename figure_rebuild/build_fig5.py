# Rebuild the two submission TikZ sources for Fig. 5 from the extracted archive.
from pathlib import Path
import sys

root = Path.cwd()
sys.path.insert(0, str(root))
from build_severe_results import current_heatmap, resilience_plot, results

rows, _ = results()
figures = {
    'severe_current_heatmap.tikz': current_heatmap(rows),
    'severe_urgent_delays.tikz': resilience_plot(rows),
}
for name, source in figures.items():
    (root / 'figures' / name).write_bytes(source.replace('\n', '\r\n').encode('utf-8'))
