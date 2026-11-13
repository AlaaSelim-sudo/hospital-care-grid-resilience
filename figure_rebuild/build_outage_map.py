# Rebuild the submission map from archived county aggregates and boundaries.
import argparse
import csv
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('study', type=Path, help='Extracted reproduction.zip directory')
args = parser.parse_args()
F = args.study / 'figures'
F.mkdir(exist_ok=True)
with Path(__file__).with_name('County_Event_Power.csv').open(encoding='utf-8-sig', newline='') as f:
    C = list(csv.DictReader(f))
G = json.loads((args.study / 'data' / 'counties_NH_VT_ME.geojson').read_text(encoding='utf-8'))['features']
def polygon(feature,xy,col,hatch=False):
 g=feature['geometry'];polys=g['coordinates'] if g['type']=='MultiPolygon' else [g['coordinates']];s=''
 for poly in polys:
  paths=[' -- '.join(f'({xy(v)[0]:.4f},{xy(v)[1]:.4f})' for v in ring)+' -- cycle' for ring in poly]
  s+=r'\path[fill='+col+r',draw=white,line width=.22pt,even odd rule] '+' '.join(paths)+';\n'
  if hatch:s+=r'\path[pattern=north east lines,pattern color=navy!28,even odd rule] '+' '.join(paths)+';\n'
 return s
def text(x,y,t,options=''):return rf'\node[font=\scriptsize{options}] at ({x},{y}) {{{t}}};'+'\n'
maps=r'\begin{tikzpicture}[x=1cm,y=1cm]\path[use as bounding box](0,0)rectangle(17.6,5.45);'+'\n'
for ox,case,title in [(0,'july2023','July 9--15, 2023'),(5.85,'dec2023','December 18--24, 2023'),(11.7,'apr2024','April 3--9, 2024')]:
 maps+=text(ox+2.65,5.2,title)
 xy=lambda v:(ox+.25+.78*.73*(v[0]+73.5),1.03+.78*(v[1]-42.65))
 for f in G:
  row=next(r for r in C if r['case']==case and r['fips']==f['properties']['GEOID']);val=row['recorded_customer_outage_hours'];k=min(4,max(0,math.floor(math.log10(max(1,float(val))))-2)) if val!='' else None
  maps+=polygon(f,xy,'unknown' if k is None else 'c'+str(k),float(row['recorded_coverage_pct'])<75 and k is not None)
 for label,coord in [('VT',(-73.0,42.49)),('NH',(-71.3,42.49)),('ME',(-68.2,43.88))]:
  x,y=xy(coord);maps+=text(x,y,label,',text=teal')
for i,label in enumerate(['$<10^3$','$10^3$--$10^4$','$10^4$--$10^5$','$10^5$--$10^6$',r'$\ge10^6$']):
 x=1.0+3.1*i;maps+=rf'\fill[c{i}]({x},.38)rectangle({x+3.1},.60);'+text(x+1.55,.17,label)
maps+=r'\end{tikzpicture}'
(F/'county_outage_maps.tikz').write_bytes(maps.replace('\n','\r\n').encode('utf-8'))
