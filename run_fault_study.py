"""Joint cohort-care, backup energy, and post-fault AC-constrained restoration.
County observations inform demand/timing; equipment/fault/repair scenarios are assumed.
"""
from pathlib import Path
import json,csv,math,time,sys
from datetime import datetime,timedelta
from collections import defaultdict
from fault_network import Network,background,gc,np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
P=Path(__file__).resolve().parent
H=72;ARRIVAL_H=48;DT=.25;KQ=gc.KQ
def read(n):return list(csv.DictReader((P/n).open(encoding='utf-8-sig')))
def write(n,rr):
 if not rr:return
 with (P/n).open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rr for k in r)));w.writeheader();w.writerows(rr)
def scenario(case,repair=12,tie_delay=6,rating=600,tie_rating=500):
 reserve=[r for r in read('data/'+case+'_EAGLEI_NH_VT_ME.csv') if r['fips_code']=='23011']
 peak=max(reserve,key=lambda r:float(r.get('customers_out',r.get('sum'))))
 event_start=datetime.fromisoformat('2023-12-18' if case=='dec2023' else '2024-04-03')
 local=datetime.fromisoformat(peak['run_start_time'])+timedelta(hours=-5 if case=='dec2023' else -4)
 peak_hour=int((local-event_start).total_seconds()/3600);start=max(0,min(96,peak_hour-12));fault=peak_hour-start
 # Cases are windows around a county peak, not a measured fault at any hospital.
 A=np.stack([np.load(P/f'clinical_{case}_H{h+1}.npz')['A'][:,start:start+ARRIVAL_H] for h in range(3)])
 stages=['normal' if t<fault or t>=fault+repair else ('island' if t<fault+tie_delay else 'backfeed') for t in range(H)]
 prior=read('Event_Hospital_Dispatch.csv');ref=np.zeros((3,H))
 for h in range(3):
  rr={int(float(r['hour'])):float(r['care_kw']) for r in prior if r['case']==case and r['policy']=='no_flex' and r['hospital']==f'H{h+1}'}
  for t in range(H):ref[h,t]=rr.get(start+t,60.) if t<ARRIVAL_H+12 else 60.
 return dict(case=case,start_hour=start,fault_hour=fault,fault_UTC=peak['run_start_time'],recorded_peak_customers_out=float(peak.get('customers_out',peak.get('sum'))),repair_h=repair,tie_delay_h=tie_delay,xfm_kva=rating,tie_rating_A=tie_rating,A=A,stages=stages,ref=ref,window_start=(event_start+timedelta(hours=start)).isoformat())
def operating_background(s,t,k):
 scale=background(s['start_hour']+t+k*.25)
 return scale*s.get('backfeed_background_fraction',1.) if s['stages'][t]=='backfeed' else scale
def linear_data(s):
 key=f"v3_{s['case']}_r{s['repair_h']}_d{s['tie_delay_h']}_x{s['xfm_kva']}_a{s['tie_rating_A']}"
 if s['repair_h']>=24:key='stress_reference_'+key
 if 'backfeed_background_fraction' in s:key+=f"_pickup{s['backfeed_background_fraction']:g}"
 cache=P/'cache'/f'{key}.npz'
 if cache.exists():return list(np.load(cache,allow_pickle=True)['data'])
 f=Network(rating=s['xfm_kva'],tie_rating=s['tie_rating_A']);data=[]
 for t in range(H):
  stage=s['stages'][t];p=s['ref'][:,t].copy();q=KQ*p
  if stage=='island':p[f.affected]=0;q[f.affected]=0
  if stage=='backfeed' and s['repair_h']>=24:p[:]=0;q[:]=-50
  if stage=='backfeed' and 'backfeed_background_fraction' in s:q[:]=-25
  for k in range(4):
   z,dv,di=f.sensitivity(stage,operating_background(s,t,k),p,q)
   data.append((t,k,z,dv,di,p.copy(),q.copy()))
 cache.parent.mkdir(exist_ok=True);np.savez_compressed(cache,data=np.array(data,dtype=object));return data
class Model:
 def __init__(self,s,policy,linear,tighten_v=0.,tighten_i=0.):
  self.s=s;self.policy=policy;self.linear=linear;self.rows=[];self.lo=[];self.hi=[];self.lb=[];self.ub=[];self.integer=[];self.names={};self.jobs=[]
  def var(name,shape,low=0,high=np.inf,integer=False):
   size=int(np.prod(shape));ix=np.arange(len(self.lb),len(self.lb)+size).reshape(shape);self.names[name]=ix;self.lb.extend(np.broadcast_to(low,shape).ravel());self.ub.extend(np.broadcast_to(high,shape).ravel());self.integer.extend([int(integer)]*size);return ix
  self.var=var
  y=var('plant',(3,H),0,1,True);load=var('load',(3,H),0,250);p=var('p',(3,H),0,300);q=var('q',(3,H),-100,100)
  c=var('charge',(3,H),0,100);d=var('discharge',(3,H),0,100);e=var('energy',(3,H+1),20,200);b=var('battery_mode',(3,H),0,1,True)
  gen=var('generator',(3,H),0,100);on=var('generator_on',(3,H),0,1,True);qs=var('q_support',(3,H),-50,50)
  self.occ=[[defaultdict(float) for t in range(H)] for h in range(3)];self.start=[[{} for t in range(H)] for h in range(3)]
  for h in range(3):
   for k,dur in enumerate([4,2]):
    for r in range(ARRIVAL_H):
     n=int(s['A'][h,k,r]);choices=[]
     if not n:continue
     for st in range(r,min(H-dur,r+12)+1):
      ix=int(var('job'+str(len(self.jobs)),(1,),0,n,True)[0]);choices.append(ix);self.jobs.append((ix,h,k,r,st,dur))
      for t in range(st,st+dur):self.occ[h][t][ix]=k
      if k==0:self.start[h][st][ix]=10.
     self.row({ix:1 for ix in choices},hi=n)
   self.lb[e[h,0]]=self.ub[e[h,0]]=200;self.lb[e[h,-1]]=160
   self.row({**{int(gen[h,t]):.30 for t in range(H)},**{int(on[h,t]):1.5 for t in range(H)}},hi=s.get('fuel_limit_l',(160.,)*3)[h])
   for t in range(H):
    fixedA=4 if t<2 else 0;fixedB=8 if t<1 else 0
    if fixedA+fixedB:self.lb[y[h,t]]=1
    terms={int(load[h,t]):1,int(y[h,t]):-60.}
    for ix,k in self.occ[h][t].items():terms[ix]=-(4 if k==0 else 1.5)
    for ix,val in self.start[h][t].items():terms[ix]=terms.get(ix,0)-val
    self.row(terms,4*fixedA+1.5*fixedB,4*fixedA+1.5*fixedB)
    self.row({**{ix:1 for ix in self.occ[h][t]},int(y[h,t]):-35},hi=-fixedA-fixedB)
    self.row({**{ix:(2 if k==0 else 1) for ix,k in self.occ[h][t].items()},int(y[h,t]):-20},hi=-2*fixedA-fixedB)
    self.row({int(p[h,t]):1,int(c[h,t]):-1,int(d[h,t]):1,int(gen[h,t]):1,int(load[h,t]):-1},0,0)
    self.row({int(q[h,t]):1,int(qs[h,t]):1,int(load[h,t]):-KQ},0,0)
    self.row({int(e[h,t+1]):1,int(e[h,t]):-1,int(c[h,t]):-.95,int(d[h,t]):1/.95},0,0)
    self.row({int(c[h,t]):1,int(b[h,t]):-100},hi=0);self.row({int(d[h,t]):1,int(b[h,t]):100},hi=100)
    self.row({int(gen[h,t]):1,int(on[h,t]):-100},hi=0);self.row({int(gen[h,t]):-1,int(on[h,t]):10},hi=0)
    for angle in (np.arange(16)+.5)*2*np.pi/16:
     co,si=np.cos(angle),np.sin(angle);self.row({int(d[h,t]):co,int(c[h,t]):-co,int(qs[h,t]):si},hi=110*np.cos(np.pi/16))
    disconnected=(s['stages'][t]=='island' or (policy=='island_only' and s['stages'][t]=='backfeed')) and h in [0,1]
    if disconnected:
     self.ub[p[h,t]]=self.lb[p[h,t]]=0;self.ub[q[h,t]]=self.lb[q[h,t]]=0
    if policy in ['autonomous','island_only'] and not disconnected:self.ub[qs[h,t]]=self.lb[qs[h,t]]=0
  if policy not in ['autonomous','island_only']:
   for t,k,z,dv,di,p0,q0 in linear:
    for j in np.where(z['energized'])[0]:
     terms={};constant=z['v'][j]
     for h in range(3):
      terms[int(p[h,t])]=float(dv[j,2*h]);terms[int(q[h,t])]=float(dv[j,2*h+1]);constant-=dv[j,2*h]*p0[h]+dv[j,2*h+1]*q0[h]
     self.row(terms,.95015+tighten_v-constant,1.04995-constant)
    if policy!='voltage_only':
     for j in range(len(z['loading'])):
      terms={};constant=z['loading'][j]
      for h in range(3):
       terms[int(p[h,t])]=float(di[j,2*h]);terms[int(q[h,t])]=float(di[j,2*h+1]);constant-=di[j,2*h]*p0[h]+di[j,2*h+1]*q0[h]
      self.row(terms,hi=.998-tighten_i-constant)
  self.n=len(self.lb);ii=[];jj=[];vv=[]
  for i,r in enumerate(self.rows):
   for j,v in r.items():
    if abs(v)>1e-12:ii.append(i);jj.append(j);vv.append(v)
  self.M=coo_matrix((vv,(ii,jj)),shape=(len(self.rows),self.n)).tocsc()
 def row(self,terms,lo=-np.inf,hi=np.inf):self.rows.append(terms);self.lo.append(lo);self.hi.append(hi)
 def optimize(self,seconds=35):
  bounds=Bounds(self.lb,self.ub);constraints=[LinearConstraint(self.M,self.lo,self.hi)];logs=[]
  critical=np.zeros(self.n);critical[self.names['plant']]=-1
  z=milp(critical,integrality=self.integer,bounds=bounds,constraints=constraints,options={'time_limit':seconds,'mip_rel_gap':1e-7})
  logs.append(self.log(z,'critical_service_hours'))
  if z.x is None:return None,logs
  # Equal-priority essential service first under every baseline.
  best=round(-float(critical@z.x));constraints.append(LinearConstraint(critical,-best-1e-7,-best+1e-7))
  objective=np.zeros(self.n)
  for ix,h,k,r,st,dur in self.jobs:
   if self.policy=='energy_first':objective[ix]=-(dur*(4 if k==0 else 1.5)+(10 if k==0 else 0))
   else:objective[ix]=-(4 if k==0 else 1)*(100*(st-r<=[2,8][k])+10-.01*(st-r))
  primary=objective.copy();cost=np.zeros(self.n);care_tie=np.zeros(self.n)
  for ix,h,k,r,st,dur in self.jobs:care_tie[ix]=-(4 if k==0 else 1)*(100*(st-r<=[2,8][k])+10-.01*(st-r))
  for h in range(3):
   for t in range(H):
    hour=(self.s['start_hour']+t)%24;price=.4 if 16<=hour<19 else (.12 if hour<7 or hour>=22 else .2)
    cost[self.names['p'][h,t]]=price;cost[self.names['generator'][h,t]]=.3*1.2;cost[self.names['generator_on'][h,t]]=1.5*1.2
    cost[self.names['charge'][h,t]]=cost[self.names['discharge'][h,t]]=.02
    objective[self.names['q_support'][h,t]]+=1e-7
  objective+=cost*1e-5
  z2=milp(objective,integrality=self.integer,bounds=bounds,constraints=constraints,options={'time_limit':seconds,'mip_rel_gap':1e-6})
  logs.append(self.log(z2,'energy' if self.policy=='energy_first' else 'care'))
  if z2.x is None:return z.x,logs
  # Tertiary cost at the obtained service objective, without claiming an unclosed gap.
  target=float(primary@z2.x);constraints.append(LinearConstraint(primary,-np.inf,target+1e-5))
  if self.policy=='energy_first':
   tie=milp(care_tie+cost*1e-5,integrality=self.integer,bounds=bounds,constraints=constraints,options={'time_limit':seconds,'mip_rel_gap':1e-6})
   logs.append(self.log(tie,'clinical_tie_break'))
   if tie.x is not None:
    z2=tie;constraints.append(LinearConstraint(care_tie,-np.inf,float(care_tie@tie.x)+1e-5))
  z3=milp(cost,integrality=self.integer,bounds=bounds,constraints=constraints,options={'time_limit':seconds,'mip_rel_gap':1e-5})
  logs.append(self.log(z3,'cost'))
  return (z3.x if z3.x is not None else z2.x),logs
 @staticmethod
 def log(z,stage):
  return dict(stage=stage,status=int(z.status),message=z.message,objective=None if z.fun is None else float(z.fun),gap=getattr(z,'mip_gap',None),dual_bound=getattr(z,'mip_dual_bound',None))
def evaluate(s,m,x):
 vals={n:x[ix] for n,ix in m.names.items() if not n.startswith('job')};f=Network(rating=s['xfm_kva'],tie_rating=s['tie_rating_A']);times=[];nodes=[];loading=[];care=[]
 totals=dict(arrivals=int(s['A'].sum()),completed=0,timely=0,high_timely=0,high_arrivals=int(s['A'][:,0].sum()),waiting_patient_h=0.,critical_unserved_kwh=float((1-vals['plant']).sum()*60))
 for ix,h,k,r,st,dur in m.jobs:
  n=int(round(x[ix]))
  if not n:continue
  totals['completed']+=n;totals['timely']+=n*(st-r<=[2,8][k]);totals['high_timely']+=n*(k==0 and st-r<=2);totals['waiting_patient_h']+=n*(st-r)
  care.append(dict(hospital=f'H{h+1}',care_class=k,arrival_hour=r,start_hour=st,duration_h=dur,patients=n))
 for t in range(H):
  for k in range(4):
   if m.policy=='island_only' and s['stages'][t]=='backfeed':
    f.topology('backfeed')
    for h in [1,2]:
     gc.dss.Text.Command(f'edit line.H{h}Supply enabled=no');gc.dss.Text.Command(f'edit load.H{h} enabled=no')
   z=f.solve_state(s['stages'][t],operating_background(s,t,k),vals['p'][:,t],vals['q'][:,t]);vv=z['v'][z['energized']];lc=z['loading'];ep=int(np.sum(vv<.95-1e-7)+np.sum(vv>1.05+1e-7))
   hp=[];risk=0.
   for h in range(3):
    connected=not((s['stages'][t]=='island' or (m.policy=='island_only' and s['stages'][t]=='backfeed')) and h in [0,1]);ids=[j for j,n in enumerate(z['nodes']) if n.split('.')[0]==f'h{h+1}'];vmin=float(min(z['v'][ids])) if connected else None;hp.append(vmin)
    # Power exposure is a separate risk screen, not a retrospective cancellation of patients.
    if connected and vmin<.95-1e-7:risk+=vals['load'][h,t]*DT
   local=(datetime.fromisoformat(s['window_start'])+timedelta(hours=t+k*.25)).isoformat()
   times.append(dict(hour=t+k*.25,local_time=local,stage=s['stages'][t],min_energized_voltage=vv.min(),max_energized_voltage=vv.max(),violating_node_phases=ep,max_current_ratio=lc.max(),overloaded_current_entries=int(sum(lc>1+1e-7)),care_kw=vals['load'][:,t].sum(),hospital_grid_kw=vals['p'][:,t].sum(),generator_kw=vals['generator'][:,t].sum(),battery_kw=(vals['discharge'][:,t]-vals['charge'][:,t]).sum(),q_support_kvar=vals['q_support'][:,t].sum(),SOC_kwh=vals['energy'][:,t].sum(),energized_service_fraction=float(vals['plant'][:,t].mean()),undervoltage_exposed_hospital_kwh=risk))
   for name,v,on in zip(z['nodes'],z['v'],z['energized']):nodes.append(dict(hour=t+k*.25,node=name,voltage_pu=v,energized=bool(on),stage=s['stages'][t]))
   groups=defaultdict(list)
   for name,v in zip(z['current_keys'],lc):groups[name.split(':')[0]].append(v)
   for name,v in groups.items():loading.append(dict(hour=t+k*.25,element=name,max_current_ratio=max(v),stage=s['stages'][t]))
 residual=max(float(np.max(np.maximum(np.asarray(m.lo)-m.M@x,0))),float(np.max(np.maximum(m.M@x-np.asarray(m.hi),0))))
 totals.update(timely_pct=100*totals['timely']/totals['arrivals'],completed_pct=100*totals['completed']/totals['arrivals'],high_timely_pct=100*totals['high_timely']/totals['high_arrivals'],service_continuity_pct=100*(1-totals['critical_unserved_kwh']/(3*H*60)),min_voltage=min(r['min_energized_voltage'] for r in times),max_voltage=max(r['max_energized_voltage'] for r in times),voltage_violation_h=sum(r['violating_node_phases']>0 for r in times)*DT,node_phase_violation_h=sum(r['violating_node_phases'] for r in times)*DT,max_current_ratio=max(r['max_current_ratio'] for r in times),overload_h=sum(r['overloaded_current_entries']>0 for r in times)*DT,exposed_hospital_kwh=sum(r['undervoltage_exposed_hospital_kwh'] for r in times),fuel_litres=float((.3*vals['generator']+1.5*vals['generator_on']).sum()),care_kwh=float(vals['load'].sum()),linear_residual=residual)
 prices=np.array([.4 if 16<=(s['start_hour']+t)%24<19 else (.12 if (s['start_hour']+t)%24<7 or (s['start_hour']+t)%24>=22 else .2) for t in range(H)])
 totals['operating_cost_dollars']=float(np.sum(vals['p']*prices)+totals['fuel_litres']*1.2+.02*np.sum(vals['charge']+vals['discharge']))
 totals['AC_voltage_thermal_feasible']=bool(totals['voltage_violation_h']==0 and totals['overload_h']==0)
 return totals,times,nodes,loading,care,vals
def main():
 policies=['island_only','autonomous','voltage_only','energy_first','coordinated'];results=[];alltimes=[];allnodes=[];allloading=[];allcare=[];logs=[];specs=[]
 seconds=180 if '--refine' in sys.argv or '--long-refine' in sys.argv else 35
 cases=[scenario(case) for case in ['dec2023','apr2024']]
 if '--sensitivity' in sys.argv:
  cases=[scenario('dec2023',repair=repair,tie_delay=delay,rating=kva) for repair,delay,kva in [(12,6,500),(12,6,650),(18,12,600),(24,18,600)]];policies=['autonomous','coordinated']
 prefix='Sensitivity_' if '--sensitivity' in sys.argv else 'Fault_'
 for s in cases:
  key=f"v4_{s['case']}_r{s['repair_h']}_d{s['tie_delay_h']}_x{s['xfm_kva']}_a{s['tie_rating_A']}";linear=linear_data(s);specs.append({k:v for k,v in s.items() if k not in ['A','ref']})
  for policy in policies:
   path=P/'solutions'/f'{key}_{policy}.npz';path.parent.mkdir(exist_ok=True)
   use_saved=path.exists() and '--resolve' not in sys.argv
   if use_saved and '--refine' in sys.argv:
    prior_logs=json.loads(path.with_suffix('.json').read_text())
    use_saved=not any(r['objective'] is None or (r.get('gap') is not None and r['gap']>.01) for r in prior_logs)
   if use_saved:
    saved=np.load(path);m=Model(s,policy,linear,float(saved['tighten_v']),float(saved['tighten_i']));x=saved['x'];ll=json.loads(path.with_suffix('.json').read_text())
   else:
    print('SOLVING',key,policy,flush=True);tv=0. if s['repair_h']>=24 else .0003;ti=.01
    for attempt in range(4):
     m=Model(s,policy,linear,tv,ti);x,ll=m.optimize(seconds)
     if x is None:print('INFEASIBLE',key,policy,ll,flush=True);break
     ev,tt,nn,ii,cc,vv=evaluate(s,m,x)
     validv=policy in ['autonomous','island_only'] or(ev['min_voltage']>=.95-1e-7 and ev['max_voltage']<=1.05+1e-7)
     validi=policy in ['autonomous','voltage_only','island_only'] or ev['max_current_ratio']<=1+1e-7
     if validv and validi:break
     if not validv:tv+=max(.0003,.95-ev['min_voltage']+.0001,ev['max_voltage']-1.05+.0001)
     if not validi:ti+=max(.002,ev['max_current_ratio']-1+.001)
     print('AC_TIGHTEN',attempt,tv,ti,flush=True)
    if x is None:logs.append(dict(scenario=key,policy=policy,solver=ll));continue
    np.savez_compressed(path,x=x,tighten_v=tv,tighten_i=ti,**vv);path.with_suffix('.json').write_text(json.dumps(ll,indent=2))
   ev,tt,nn,ii,cc,vv=evaluate(s,m,x);ev.update(scenario=key,case=s['case'],policy=policy,repair_h=s['repair_h'],island_h=min(s['tie_delay_h'],s['repair_h']),tie_rating_A=s['tie_rating_A'],xfm_kva=s['xfm_kva']);results.append(ev)
   for rows,out in [(tt,alltimes),(nn,allnodes),(ii,allloading),(cc,allcare)]:out.extend(dict(scenario=key,case=s['case'],policy=policy,**r) for r in rows)
   logs.append(dict(scenario=key,policy=policy,solver=ll));print('RESULT',json.dumps(ev),flush=True)
   for name,rr in [('Results',results),('Time_Series',alltimes),('Node_Voltages',allnodes),('Current_Loading',allloading),('Treatment_Starts',allcare)]:write(prefix+name+'.csv',rr)
   (P/(prefix+'Solver_Log.json')).write_text(json.dumps(logs,indent=2))
 (P/(prefix+'Design.json')).write_text(json.dumps(specs,indent=2))
if __name__=='__main__':main()
