"""Native TikZ figures with continuous, proportional colorbars and fixed scales."""
from run_fault_study import read,write,P,scenario,np,H
from collections import defaultdict
import math,json
F=P/'figures'
VST=[(.90,'AB581D'),(.95,'F4E8D9'),(1.,'69A7B2'),(1.05,'193650')]
IST=[(0,'FFFFFF'),(.8,'A4CDD0'),(1.,'F4E8D9'),(1.2,'AB581D')]
DST=[(-1,'007F86'),(0,'FFFFFF'),(1,'C57939')]
BUS=['sourcebus','650','rg60','632','645','646','633','634','670','671','680','684','611','652','692','675','h1','h2','h3']
ELEMENTS=['Line.650632','Line.632670','Line.670671','Line.671680','Line.632633','Line.632645','Line.645646','Line.692675','Line.671684','Line.684611','Line.684652','Line.671692','Line.h1supply','Line.h2supply','Line.h3supply','Line.restorationtie','Transformer.sub','Transformer.reg1','Transformer.reg2','Transformer.reg3','Transformer.xfm1']
def color(v,stops):
 v=max(stops[0][0],min(stops[-1][0],v))
 for (a,ca),(b,cb) in zip(stops,stops[1:]):
  if a<=v<=b:
   u=(v-a)/(b-a);return ''.join(f'{round((1-u)*int(ca[k:k+2],16)+u*int(cb[k:k+2],16)):02X}' for k in [0,2,4])
 return stops[-1][1]
def rect(x,y,w,h,c):return rf'\definecolor{{cell}}{{HTML}}{{{c}}}\fill[cell]({x},{y})rectangle({x+w},{y+h});'
def colorbar(x,y,width,stops,ticks,label,fmt=lambda x:str(x)):
 a,b=stops[0][0],stops[-1][0];s=''
 for k in range(120):s+=rect(x+k*width/120,y,width/120+.001,.65,color(a+(b-a)*(k+.5)/120,stops))
 s+=rf'\draw[gray]({x},{y})rectangle({x+width},{y+.65});'
 for v in ticks:
  xx=x+width*(v-a)/(b-a);s+=rf'\draw({xx},{y})--({xx},{y-.22});\node[anchor=north]at({xx},{y-.3}){{{fmt(v)}}};'
 threshold=.95 if a==.90 else (1. if a==0 and b==1.2 else None)
 if threshold is not None:
  xx=x+width*(threshold-a)/(b-a);s+=rf'\draw[red,thick]({xx},{y-.1})--({xx},{y+.75});'
 s+=rf'\node[anchor=north]at({x+width/2},{y-1.5}){{{label}}};';return s
def label(n):return {'sourcebus':'Source','rg60':'REG','h1':'H1','h2':'H2','h3':'H3','Transformer.sub':'Sub','Transformer.reg1':'Reg A','Transformer.reg2':'Reg B','Transformer.reg3':'Reg C','Transformer.xfm1':'XFM1','Line.restorationtie':'Tie'}.get(n,n.replace('Line.','').replace('supply',''))
def axis(x,y,nrows,s):
 r=rf'\draw[gray]({x},{y})rectangle({x+72},{y+nrows});'
 for h in [0,12,24,36,48,60,72]:r+=rf'\node[anchor=north]at({x+h},{y-.2}){{{h}}};'
 r+=rf'\node[anchor=north]at({x+36},{y-1.4}){{Elapsed hour}};'
 for h in [s['fault_hour'],s['fault_hour']+s['tie_delay_h'],s['fault_hour']+s['repair_h']]:r+=rf'\draw[densely dashed,black!65]({x+h},{y})--({x+h},{y+nrows});'
 return r
def main():
 nodes=read('Fault_Node_Voltages.csv');curr=read('Fault_Current_Loading.csv');times=read('Fault_Time_Series.csv');res=read('Fault_Results.csv')
 v=defaultdict(list);ii=defaultdict(list);vq=defaultdict(list);iq={}
 for r in nodes:
  if r['energized']=='True':
   v[r['case'],r['policy'],r['node'].split('.')[0],int(float(r['hour']))].append(float(r['voltage_pu']))
   vq[r['case'],r['policy'],r['node'].split('.')[0],float(r['hour'])].append(float(r['voltage_pu']))
 for r in curr:
  ii[r['case'],r['policy'],r['element'],int(float(r['hour']))].append(float(r['max_current_ratio']))
  iq[r['case'],r['policy'],r['element'],float(r['hour'])]=float(r['max_current_ratio'])
 hourly=[];delta_rows=[]
 for case in ['dec2023','apr2024']:
  s=scenario(case)
  tex=r'\begin{tikzpicture}[x=.08cm,y=.24cm,every node/.style={font=\fontsize{6}{6.5}\selectfont}]'
  for col,pol in enumerate(['autonomous','coordinated']):
   x=col*100;y=33
   tex+=rf'\node[font=\scriptsize]at({x+36},{y+21.5}){{'+('Local dispatch' if col==0 else 'Joint care--utility coordination')+'};'
   for j,b in enumerate(BUS):
    yy=y+18-j;tex+=rf'\node[anchor=east]at({x-1},{yy+.5}){{{label(b)}}};'
    for h in range(72):
     arr=v[case,pol,b,h];val=min(arr) if arr else None;c='DADDE1' if val is None else color(val,VST);tex+=rect(x+h,yy,1,1,c)
     hourly.append(dict(case=case,policy=pol,position=b,kind='minimum_voltage_pu',hour=h,value='' if val is None else val))
   tex+=axis(x,y,19,s)
   for j,b in enumerate(ELEMENTS):
    yy=20-j;tex+=rf'\node[anchor=east]at({x-1},{yy+.5}){{{label(b)}}};'
    for h in range(72):
     arr=ii[case,pol,b,h];val=max(arr) if arr else None;tex+=rect(x+h,yy,1,1,'DADDE1' if val is None else color(val,IST))
     hourly.append(dict(case=case,policy=pol,position=b,kind='maximum_current_ratio',hour=h,value='' if val is None else val))
   tex+=axis(x,0,21,s)
  tex+=colorbar(181,29,0,VST,[], '') if False else ''
  # Both panels in each row share a continuous proportional colorbar.
  tex+=colorbar(12,28.2,144,VST,[.90,.95,1,1.05],'Minimum energized phase voltage (p.u.)',lambda x:f'{x:.2f}')
  tex+=colorbar(12,-4.5,144,IST,[0,.4,.8,1,1.2],r'Maximum phase current / scenario rating (\%)',lambda x:f'{100*x:.0f}')
  tex+=r'\node[anchor=west]at(-1,-8){Gray = isolated or open; dashed = fault / tie / repair.};\end{tikzpicture}'
  (F/f'fault_heatmap_{case}.tikz').write_text(tex.replace(';',';\n'))
  # Main paper: current maps only; the larger four-panel maps stay in supplement.
  ct=r'\begin{tikzpicture}[x=.08cm,y=.24cm,every node/.style={font=\fontsize{6}{6.5}\selectfont}]'
  for col,pol in enumerate(['autonomous','coordinated']):
   xx=col*100;ct+=rf'\node[font=\scriptsize]at({xx+36},23){{'+('Local dispatch' if col==0 else 'Joint care--utility coordination')+'};'
   for j,elem in enumerate(ELEMENTS):
    yy=20-j;ct+=rf'\node[anchor=east]at({xx-1},{yy+.5}){{{label(elem)}}};'
    for h in range(72):
     arr=ii[case,pol,elem,h];ct+=rect(xx+h,yy,1,1,'DADDE1' if not arr else color(max(arr),IST))
   ct+=axis(xx,0,21,s)
  ct+=colorbar(12,-4.5,144,IST,[0,.4,.8,1,1.2],r'Maximum phase current / scenario rating (\%)',lambda x:f'{100*x:.0f}')
  ct+=r'\end{tikzpicture}';(F/f'current_heatmap_{case}.tikz').write_text(ct.replace(';',';\n'))
  # Actual differences use a single shared scale across all rows and both cases.
  tex=r'\begin{tikzpicture}[x=.08cm,y=.2cm,every node/.style={font=\scriptsize}]'
  for col,(kind,positions,src,scale) in enumerate([('Voltage change',BUS,v,.025),('Current-loading change',ELEMENTS,ii,.40)]):
   x=col*105;n=len(positions);tex+=rf'\node at({x+36},{n+1.5}){{{kind}}};'
   for j,b in enumerate(positions):
    yy=n-1-j;tex+=rf'\node[anchor=east]at({x-1},{yy+.5}){{{label(b)}}};'
    for h in range(72):
     delta=[]
     for k in range(4):
      keya=case,'autonomous',b,h+k*.25;keyb=case,'coordinated',b,h+k*.25
      if col==0 and vq[keya] and vq[keyb]:delta.append(min(vq[keyb])-min(vq[keya]))
      if col==1 and keya in iq and keyb in iq:delta.append(iq[keyb]-iq[keya])
     val=sum(delta)/len(delta) if delta else None
     if val is not None:assert abs(val)<=scale,(case,b,h,val,scale)
     delta_rows.append(dict(case=case,position=b,kind='paired_voltage_change_pu' if col==0 else 'paired_current_ratio_change',hour=h,value='' if val is None else val,colorbar_abs_limit=scale))
     tex+=rect(x+h,yy,1,1,'DADDE1' if val is None else color(val/scale,DST))
   tex+=axis(x,0,n,s)
   tex+=colorbar(x+3,-4.5,66,[(a*scale,c) for a,c in DST],[-scale,0,scale],'After minus before (p.u.)' if col==0 else 'After minus before (percentage points)',lambda x:f'{x:.3f}' if col==0 else f'{x*100:+.0f}')
  tex+=r'\end{tikzpicture}';(F/f'fault_difference_{case}.tikz').write_text(tex.replace(';',';\n'))
 write('Plotted_Hourly_Heatmap_Cells.csv',hourly)
 write('Plotted_Hourly_Changes.csv',delta_rows)
 # One topology only. Positions follow native IEEE branches, with three added sites.
 coords={'sourcebus':(0,7),'650':(0,6.4),'rg60':(0,5.7),'632':(0,4.8),'633':(3,4.8),'634':(5,4.8),'645':(-2,4.8),'646':(-4,4.8),'670':(0,3.5),'671':(0,2.2),'692':(2,2.2),'675':(4,2.2),'684':(-2,2.2),'611':(-4,2.2),'652':(-2,.4),'680':(0,.4),'h1':(6,2.2),'h2':(0,-1.1),'h3':(7,4.8)}
 edges=[('sourcebus','650'),('650','rg60'),('rg60','632'),('632','633'),('633','634'),('632','645'),('645','646'),('632','670'),('670','671'),('671','680'),('671','692'),('692','675'),('671','684'),('684','611'),('684','652'),('675','h1'),('680','h2'),('634','h3')]
 tex=r'\begin{tikzpicture}[x=.8cm,y=.65cm,>=Latex,every node/.style={font=\scriptsize}]'
 for a,b in edges:
  sty='orange,thick' if b.startswith('h') else ('red,thick' if (a,b)==('670','671') else 'navy')
  tex+=rf'\draw[{sty}]{coords[a]}--{coords[b]};'
 for n,(x,y) in coords.items():
  if n.startswith('h'):tex+=rf'\node[draw=orange,fill=white,inner sep=3pt]at({x},{y}){{{n.upper()}}};'
  else:tex+=rf'\fill[navy]({x},{y})circle(1.3pt);\node[anchor=south west,fill=white,inner sep=.5pt]at({x+.04},{y+.05}){{{label(n)}}};'
 tex+=r'\draw[teal,dashed,thick,->](3,4.8)to[out=-45,in=0](0,.4);\node[teal,fill=white,align=center]at(3.8,.8){Normally-open tie\\633--680; 500 A};'
 tex+=r'\node[red,font=\large]at(0,2.9){$\times$};\node[red,anchor=east]at(-.15,2.9){Fault / trip};'
 tex+=r'\node[fill=white]at(4,4.45){XFM1: 600 kVA};\node[align=left,anchor=west]at(1.6,6.6){Each modeled ED:\\100-kW generator + 160 L\\200-kWh / 100-kW battery\\$\pm50$ kvar support};'
 tex+=r'\end{tikzpicture}';(F/'ieee13_single.tikz').write_text(tex.replace(';',';\n'))
 # Common physical time plots: all phases, not selected best/worst nodes.
 styles=[('autonomous','orange','Local dispatch'),('voltage_only','purple','Voltage only'),('coordinated','teal','Joint V/I')]
 tex=r'\begin{tikzpicture}\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.7cm},width=.43\textwidth,height=4.2cm,xmin=10,xmax=28,xtick={12,18,24},grid=major,grid style={gray!15},tick label style={font=\scriptsize},label style={font=\scriptsize},legend style={font=\scriptsize,draw=none,legend columns=3,at={(1.15,1.16)},anchor=south}]'
 for col in [0,1]:
  tex+=r'\nextgroupplot[xlabel={Elapsed hour},ylabel={'+('Lowest energized voltage (p.u.)' if col==0 else r'Highest phase-current ratio (\%)')+'}]'
  tex+=r'\path[fill=gray!12,draw=none](rel axis cs:.1111,0)rectangle(rel axis cs:.4444,1);\node[font=\scriptsize,fill=white,inner sep=1pt]at(rel axis cs:.2777,.92){Isolated};'
  for policy,c,name in styles:
   rr=[r for r in times if r['case']=='dec2023' and r['policy']==policy and 10<=float(r['hour'])<=28]
   tex+=r'\addplot[const plot,thick,'+c+']coordinates{'+''.join(f'({r["hour"]},{float(r["min_energized_voltage"] if col==0 else r["max_current_ratio"])*(1 if col==0 else 100)})' for r in rr)+'};'
   if col==0:tex+=r'\addlegendentry{'+name+'}'
  threshold=.95 if col==0 else 100;tex+=rf'\addplot[black,dashed,forget plot]coordinates{{(10,{threshold})(28,{threshold})}};'
 tex+=r'\end{groupplot}\end{tikzpicture}';(F/'restoration_time.tikz').write_text(tex)
 # Actual saved resource actions and clinical service, rather than a second feeder.
 ct=read('Fault_Treatment_Starts.csv');resource_rows=[]
 tex=r'\begin{tikzpicture}\begin{groupplot}[group style={group size=2 by 2,horizontal sep=1.8cm,vertical sep=2cm},width=.43\textwidth,height=4.5cm,grid=major,grid style={gray!15},tick label style={font=\scriptsize},label style={font=\scriptsize},title style={font=\small},legend style={font=\scriptsize,draw=none,fill=white,at={(.5,.98)},anchor=north,legend columns=1}]'
 saved={p:np.load(P/'solutions'/f'v4_dec2023_r12_d6_x600_a500_{p}.npz') for p in ['autonomous','coordinated','island_only']}
 for pol,z in saved.items():
  used=np.cumsum(.3*z['generator']+1.5*z['generator_on'],axis=1)
  for h in range(3):
   for t in range(72):resource_rows.append(dict(case='dec2023',policy=pol,hospital=f'H{h+1}',hour=t,care_kw=z['load'][h,t],grid_kw=z['p'][h,t],generator_kw=z['generator'][h,t],battery_kw=z['discharge'][h,t]-z['charge'][h,t],q_support_kvar=z['q_support'][h,t],SOC_kwh=z['energy'][h,t],remaining_fuel_L=160-used[h,t],essential_available=z['plant'][h,t]))
 for panel in range(4):
  tex+=r'\nextgroupplot[xlabel={Elapsed hour},'+('xmin=10,xmax=28,xtick={12,18,24},' if panel<3 else 'xmin=0,xmax=60,')+['ymin=0,ymax=550,','ymin=-350,ymax=450,','ymin=-30,ymax=240,','ymin=0,ymax=800,'][panel]+'ylabel={'+['Grid import (kW)','Onsite active power (kW)','Inverter support (kvar)','Cumulative timely starts'][panel]+'},title={'+['(a) Shared care, lower utility import','(b) Coordinated onsite action','(c) Reactive support is explicit','(d) Timely care after reconnection'][panel]+'}]'
  if panel==0:lines=[('autonomous','p','orange','Local grid import'),('coordinated','p','teal','Joint grid import'),('coordinated','load','gray','Scheduled care demand')]
  elif panel==1:lines=[('coordinated','generator','teal','Generators'),('coordinated','battery','orange','Battery net output')]
  elif panel==2:lines=[('autonomous','q_support','orange','Local'),('coordinated','q_support','teal','Joint')]
  else:lines=[('island_only','patients','orange','Island hold'),('coordinated','patients','teal','Joint')]
  for pol,kind,col,name in lines:
   z=saved[pol]
   if kind=='patients':
    values=[sum(int(r['patients']) for r in ct if r['case']=='dec2023' and r['policy']==pol and int(r['start_hour'])<=t and int(r['start_hour'])-int(r['arrival_hour'])<=[2,8][int(r['care_class'])]) for t in range(72)]
   elif kind=='battery':values=(z['discharge']-z['charge']).sum(axis=0)
   else:values=z[kind].sum(axis=0)
   tex+=r'\addplot[const plot,thick,'+col+']coordinates{'+''.join(f'({t},{float(val)})' for t,val in enumerate(values))+'};'+r'\addlegendentry{'+name+'}'
 tex+=r'\end{groupplot}\end{tikzpicture}';(F/'hospital_utility_actions.tikz').write_text(tex)
 write('Plotted_Resource_Dispatch.csv',resource_rows)
 print('Generated one topology, continuous-colorbar heatmaps, paired changes and time traces.')
if __name__=='__main__':main()
