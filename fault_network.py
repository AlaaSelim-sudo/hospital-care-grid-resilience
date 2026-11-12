"""Assumed fault/isolation/backfeed topology on the preserved IEEE-13 source."""
import grid_core as gc
import numpy as np,math
dss=gc.dss
LINE_RATINGS={'650632':700,'632670':600,'670671':600,'671680':600,'632633':600,'632645':300,'645646':300,'692675':400,'671684':300,'684611':300,'684652':300,'671692':600,'h1supply':120,'h2supply':120,'h3supply':300,'restorationtie':500}
class Network(gc.Feeder):
 def __init__(self,kind='trunk',rating=600,tie_rating=500):
  super().__init__();self.kind=kind;self.rating=rating;self.ratings=dict(LINE_RATINGS,restorationtie=tie_rating)
  for name,tap in zip(['reg1','reg2','reg3'],[1.0495,1.0375,1.0495]):
   dss.Transformers.Name(name);dss.Transformers.Wdg(2);dss.Transformers.Tap(tap)
  dss.Text.Command(f'edit transformer.XFM1 kvas=[{rating} {rating}]')
  # Shared benchmark normalization: retain total 634 P/Q, balance its phases.
  # Prevent the original native phase-A overload from confounding surge effects.
  for phase in ['a','b','c']:dss.Text.Command(f'edit load.634{phase} kw={400/3} kvar={290/3}')
  self.base=[(n,400/3,290/3) if n.lower() in ['634a','634b','634c'] else (n,p,q) for n,p,q in self.base]
  dss.Text.Command('edit line.H3Supply bus1=634.1.2.3 length=.005')
  dss.Text.Command('edit load.H3 kv=.48')
  dss.Text.Command('new line.RestorationTie phases=3 bus1=633.1.2.3 bus2=680.1.2.3 linecode=mtx601 length=.05 units=mi enabled=no')
  dss.Text.Command('calcvoltagebases');self.stage=None
 def topology(self,stage):
  if self.stage==stage:return
  self.stage=stage
  damaged='670671' if self.kind=='trunk' else '632633'
  dss.Text.Command('edit line.'+damaged+' enabled='+('yes' if stage=='normal' else 'no'))
  dss.Text.Command('edit line.RestorationTie enabled='+('yes' if stage=='backfeed' else 'no'))
  off=['671','692','675a','675b','675c','611','652'] if self.kind=='trunk' else ['634a','634b','634c']
  for name in off:dss.Text.Command('edit load.'+name+' enabled='+('no' if stage=='island' else 'yes'))
  if self.kind=='trunk':
   for name in ['Cap1','Cap2']:dss.Text.Command('edit capacitor.'+name+' enabled='+('no' if stage=='island' else 'yes'))
  for h in range(3):dss.Text.Command(f'edit load.H{h+1} enabled='+('no' if stage=='island' and h in self.affected else 'yes'))
  for h in range(3):dss.Text.Command(f'edit line.H{h+1}Supply enabled=yes')
  # Recreated PCC buses need their bases reset after open/close operations.
  # Otherwise raw volts can be divided by an uninitialized base after repair.
  dss.Text.Command('calcvoltagebases')
 @property
 def affected(self):return [0,1] if self.kind=='trunk' else [2]
 def solve_state(self,stage,scale,p,q):
  self.topology(stage);z=super().solve(scale,p,q)
  # Source reachability, not residual voltages on floating disconnected buses,
  # defines energized grid locations. Internal hospital islands are separate.
  edges=[]
  for cls,names in [('Line',dss.Lines.AllNames()),('Transformer',dss.Transformers.AllNames())]:
   for name in names:
    dss.Circuit.SetActiveElement(cls+'.'+name)
    if dss.CktElement.Enabled():
     buses=[b.split('.')[0].lower() for b in dss.CktElement.BusNames()]
     if len(buses)>1:edges.append((buses[0],buses[1]))
  reach={'sourcebus'}
  for _ in range(30):
   before=len(reach)
   for a,b in edges:
    if a in reach or b in reach:reach.update([a,b])
   if len(reach)==before:break
  v=np.array(z['v']);nodes=z['names'];values=[];keys=[];amps=[]
  for name in dss.Lines.AllNames():
   dss.Circuit.SetActiveElement('Line.'+name)
   if not dss.CktElement.Enabled():continue
   if not all(b.split('.')[0].lower() in reach for b in dss.CktElement.BusNames()):continue
   aa=np.array(dss.CktElement.CurrentsMagAng()[::2]);nn=dss.CktElement.NodeOrder()
   for j,(a,n) in enumerate(zip(aa,nn)):
    if n:values.append(a/self.ratings[name]);keys.append('Line.'+name+':'+str(j));amps.append(a)
  for name in dss.Transformers.AllNames():
   dss.Transformers.Name(name);ph=dss.CktElement.NumPhases();nc=dss.CktElement.NumConductors();aa=np.array(dss.CktElement.CurrentsMagAng()[::2]);nn=dss.CktElement.NodeOrder()
   if not all(b.split('.')[0].lower() in reach for b in dss.CktElement.BusNames()):continue
   for w in range(dss.Transformers.NumWindings()):
    dss.Transformers.Wdg(w+1);base=dss.Transformers.kVA()/(dss.Transformers.kV()*(math.sqrt(3) if ph==3 else 1))
    for j in range(nc):
     k=w*nc+j
     if nn[k]:values.append(aa[k]/base);keys.append('Transformer.'+name+':'+str(k));amps.append(aa[k])
  energized=np.array([n.split('.')[0].lower() in reach for n in nodes])
  return dict(v=v,nodes=nodes,loading=np.array(values),current_keys=keys,amps=np.array(amps),energized=energized,loss_kw=z['loss_kw'])
 def sensitivity(self,stage,scale,p,q):
  z=self.solve_state(stage,scale,p,q);dv=[];di=[]
  for h in range(3):
   for which in ['P','Q']:
    pp=p.copy();qq=q.copy()
    if which=='P':pp[h]+=1
    else:qq[h]+=1
    r=self.solve_state(stage,scale,pp,qq);dv.append(r['v']-z['v']);di.append(r['loading']-z['loading'])
  return z,np.array(dv).T,np.array(di).T

def background(hour):
 day=np.array([.96,.95,.94,.94,.95,.96,.97,.98,.99,1,1,1,1,1,1,1.01,1.02,1.02,1.02,1.01,1,.99,.98,.97])-.01
 return float(np.interp(hour%24,np.arange(25),np.r_[day,day[0]]))
