"""Cohort ED scheduling and modular diesel/battery MILP in physical units.
Conditional planning experiment, NOT a calibrated hospital digital twin.
"""
from pathlib import Path
import sys,json,time
import numpy as np
try:from scipy.optimize import milp, Bounds, LinearConstraint
except ImportError:
 sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'research/.deps'))
 from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix,vstack

class EDModel:
 def __init__(self,arrivals,available,cfg):
  self.A=np.asarray(arrivals,int);self.a=np.asarray(available,float);self.cfg=cfg;self.T=len(available);T=self.T
  self.lb=[];self.ub=[];self.integ=[];self.idx={};self.rows=[];self.rlo=[];self.rhi=[]
  def var(name,size,upper,integer=False):
   z=np.arange(len(self.lb),len(self.lb)+size);self.lb.extend([0.]*size);self.ub.extend(np.broadcast_to(upper,(size,)));self.integ.extend([int(integer)]*size);self.idx[name]=z;return z
  def row(terms,lo=-np.inf,hi=np.inf):
   self.rows.append(dict(terms));self.rlo.append(lo);self.rhi.append(hi)
  self.row=row;self.jobs=[];self.occ=[[[] for t in range(T)] for k in range(2)];self.startA=[[] for t in range(T)]
  self.bycohort={}
  # Class 0: higher-acuity design class, duration 4 h; class 1: routine, 2 h.
  # Delays 2/8 h are experimental targets, not clinical standards.
  for k in range(2):
   dur=cfg['duration_h'][k]
   for r in range(self.A.shape[1]):
    if not self.A[k,r]:continue
    choices=[]
    for s in range(r,min(r+cfg['maximum_wait_h'],T-dur)+1):
     ix=int(var('x'+str(len(self.jobs)),1,int(self.A[k,r]),True)[0]);choices.append(ix)
     self.jobs.append((ix,k,r,s,dur));
     for t in range(s,s+dur):self.occ[k][t].append(ix)
     if k==0:self.startA[s].append(ix)
    self.bycohort[k,r]=choices;row({i:1 for i in choices},hi=int(self.A[k,r]))
  self.fixed=np.zeros((2,T))
  for k in range(2):self.fixed[k,:cfg['initial_remaining_h'][k]]=cfg['initial_patients'][k]
  N=cfg['generator_modules'];G=cfg['module_prime_kw'];B=cfg['battery_kwh'];BP=cfg['battery_kw']
  plant=var('plant',T,1,True);grid=var('grid',T,cfg['grid_kw']*self.a)
  on=var('generator_on',N*T,1,True).reshape(N,T);power=var('generator_kw',N*T,G).reshape(N,T)
  start=var('generator_start',N*T,1,True).reshape(N,T)
  ch=var('charge_kw',T,BP);dis=var('discharge_kw',T,BP);soc=var('soc_kwh',T+1,B);mode=var('battery_charge_mode',T,1,True);spill=var('spill_kw',T,N*G)
  fuel=var('fuel_litres',T+1,cfg['tank_litres'])
  self.lb[soc[0]]=self.ub[soc[0]]=B;self.lb[fuel[0]]=self.ub[fuel[0]]=cfg['initial_fuel_litres']
  self.lb[soc[-1]]=cfg['terminal_soc_fraction']*B
  self.energy_terms={};self.fuel_terms={}
  for t in range(T):
   active={i:1 for k in range(2) for i in self.occ[k][t]};fixed=sum(self.fixed[:,t])
   row(active,hi=cfg['ED_bays']-fixed)
   row({i:cfg['nurse_units'][k] for k in range(2) for i in self.occ[k][t]},hi=cfg['staff_units']-sum(self.fixed[k,t]*cfg['nurse_units'][k] for k in range(2)))
   row({**active,int(plant[t]):-cfg['ED_bays']},hi=-fixed)
   # No powered empty department: energy baseline cannot earn idle-plant energy.
   row({int(plant[t]):1,**{i:-1 for i in active}},hi=fixed)
   load={int(plant[t]):-cfg['shared_kw']}
   for k in range(2):
    for i in self.occ[k][t]:load[i]=-cfg['occupied_bay_kw'][k]
   for i in self.startA[t]:load[i]-=cfg['diagnostic_kwh_per_class0_start']
   fixed_kw=sum(self.fixed[k,t]*cfg['occupied_bay_kw'][k] for k in range(2))
   terms={**load,int(grid[t]):1,int(ch[t]):-1,int(dis[t]):1,int(spill[t]):-1}
   for n in range(N):terms[int(power[n,t])]=1
   row(terms,fixed_kw,fixed_kw)
   row({int(spill[t]):1,**{int(power[n,t]):-1 for n in range(N)}},hi=0)
   row({int(ch[t]):1,int(mode[t]):-BP},hi=0);row({int(dis[t]):1,int(mode[t]):BP},hi=BP)
   row({int(soc[t+1]):1,int(soc[t]):-1,int(ch[t]):-cfg['charge_efficiency'],int(dis[t]):1/cfg['discharge_efficiency']},0,0)
   ft={int(fuel[t+1]):1,int(fuel[t]):-1}
   for n in range(N):
    p,y,z=map(int,(power[n,t],on[n,t],start[n,t]))
    row({p:1,y:-G},hi=0);row({p:1,y:-.75*G},lo=0)
    # Exact start indicator, generators initially off.
    terms={z:1,y:-1}
    if t:terms[int(on[n,t-1])]=1
    row(terms,lo=0);row({z:1,y:-1},hi=0)
    if t:row({z:1,int(on[n,t-1]):1},hi=1)
    for j in range(t,min(T,t+cfg['minimum_up_h'])):row({int(on[n,j]):1,z:-1},lo=0)
    ft[p]=cfg['fuel_slope_litre_per_kwh'];ft[y]=cfg['fuel_intercept_litre_per_h'];ft[z]=cfg['start_fuel_litres']
    self.fuel_terms[p]=cfg['fuel_slope_litre_per_kwh'];self.fuel_terms[y]=cfg['fuel_intercept_litre_per_h'];self.fuel_terms[z]=cfg['start_fuel_litres']
   delivery=cfg.get('delivery_litres_by_hour',{}).get(str(t),0)
   row(ft,delivery,delivery)
  self.n=len(self.lb);rr=[];cc=[];vv=[]
  for j,terms in enumerate(self.rows):
   for i,v in terms.items():rr.append(j);cc.append(i);vv.append(v)
  self.M=coo_matrix((vv,(rr,cc)),shape=(len(self.rows),self.n)).tocsc();self.lower=np.array(self.rlo);self.upper=np.array(self.rhi)
  self.scores={k:np.zeros(self.n) for k in ['timely','completed','energy','burden','fuel']}
  # Energy counts patient-associated electrical energy, excludes shared/idle plant and spills.
  for ix,k,r,s,dur in self.jobs:
   self.scores['timely'][ix]=cfg['urgency_weights'][k]*(s-r<=cfg['deadline_h'][k])
   self.scores['completed'][ix]=1
   self.scores['energy'][ix]=dur*cfg['occupied_bay_kw'][k]+(cfg['diagnostic_kwh_per_class0_start'] if k==0 else 0)
   self.scores['burden'][ix]=s-r-(T-r) # constant total no-start burden omitted
  for i,v in self.fuel_terms.items():self.scores['fuel'][i]=v
  self.scores['tie']=-self.scores['fuel']-.01*self.scores['burden']

 def solve(self,policy,label,folder,limit=12):
  # Equal physical feasible set and perfect information in all policies.
  priorities={'care':['timely','completed','tie'],'energy':['energy','timely','tie'],'throughput':['completed','timely','tie']}[policy]
  mat=self.M;lo=self.lower.copy();hi=self.upper.copy();logs=[]
  # A constructive feasible no-new-start schedule protects initial active cases
  # from the full battery, then recharges on available grid. This is a witness,
  # not a solver warm start (SciPy's milp API has no x0 argument).
  x=np.zeros(self.n);soc=self.cfg['battery_kwh'];fuel=self.cfg['initial_fuel_litres'];x[self.idx['soc_kwh'][0]]=soc;x[self.idx['fuel_litres'][0]]=fuel
  for t in range(self.T):
   load=float(np.dot(self.fixed[:,t],self.cfg['occupied_bay_kw']));active=sum(self.fixed[:,t])>0
   x[self.idx['plant'][t]]=int(active);load+=self.cfg['shared_kw']*active
   if self.a[t]:
    charge=min(self.cfg['battery_kw'],(self.cfg['battery_kwh']-soc)/self.cfg['charge_efficiency']);x[self.idx['charge_kw'][t]]=charge;x[self.idx['battery_charge_mode'][t]]=1;x[self.idx['grid'][t]]=load+charge;soc+=charge*self.cfg['charge_efficiency']
   else:x[self.idx['discharge_kw'][t]]=load;soc-=load/self.cfg['discharge_efficiency']
   x[self.idx['soc_kwh'][t+1]]=soc;x[self.idx['fuel_litres'][t+1]]=fuel
  assert max(self.M@x-self.upper)<1e-5 and max(self.lower-self.M@x)<1e-5 and min(x-np.array(self.lb))>-1e-5
  # Reuse feasible witnesses from other priority policies on identical inputs.
  for f in folder.glob(label.rsplit('_',1)[0]+'_*.npz'):
   z=np.load(f);xx=z['x']
   if len(xx)==self.n and np.array_equal(z['A'],self.A) and np.array_equal(z['a'],self.a) and json.loads(str(z['cfg']))==self.cfg and self.scores[priorities[0]]@xx>self.scores[priorities[0]]@x:x=xx.copy()
  # Resource upgrades must retain a feasible base-design witness. A larger
  # resource set must not appear worse just because a solver reached its limit.
  for design in ['fuel_x2','battery_x2','extra_module']:
   if '_'+design+'_' not in label:continue
   f=folder/(label.replace('_'+design+'_','_base_').rsplit('_',1)[0]+'_care.npz')
   if not f.exists():continue
   z=np.load(f);oldcfg=json.loads(str(z['cfg']));old=EDModel(z['A'],z['a'],oldcfg);xx=np.zeros(self.n)
   for key,ids in self.idx.items():
    source=z['x'][old.idx[key]]
    if key in ['generator_on','generator_kw','generator_start']:
     xx[ids.reshape(-1,self.T)[:oldcfg['generator_modules']].flatten()]=source
    else:xx[ids]=source
   xx[self.idx['soc_kwh']]+=self.cfg['battery_kwh']-oldcfg['battery_kwh']
   xx[self.idx['fuel_litres']]+=self.cfg['initial_fuel_litres']-oldcfg['initial_fuel_litres']
   assert max(self.M@xx-self.upper)<2e-5 and max(self.lower-self.M@xx)<2e-5 and max(np.array(self.lb)-xx)<2e-5 and max(xx-np.array(self.ub))<2e-5
   if self.scores[priorities[0]]@xx>self.scores[priorities[0]]@x:x=xx
  upper={'timely':float(np.dot(self.A.sum(axis=1),self.cfg['urgency_weights'])),'completed':float(self.A.sum()),'energy':float(sum(self.A[k].sum()*(self.cfg['duration_h'][k]*self.cfg['occupied_bay_kw'][k]+(self.cfg['diagnostic_kwh_per_class0_start'] if k==0 else 0)) for k in range(2)))}
  for stage,key in enumerate(priorities):
   obj=self.scores[key];now=time.perf_counter()
   if key in upper and abs(obj@x-upper[key])<1e-5:
    value=float(obj@x);record=dict(label=label,stage=stage+1,objective=key,status=0,value=value,gap=0.,seconds=0.,dual_bound=-upper[key],message='Feasible witness attains independent absolute domain upper bound',certificate_kind='analytic_domain_bound',used_previous_feasible_incumbent=True)
   else:
    z=milp(-obj,integrality=np.array(self.integ),bounds=Bounds(self.lb,self.ub),constraints=LinearConstraint(mat,lo,hi),options={'time_limit':limit,'mip_rel_gap':1e-5})
    old=float(obj@x);used=z.x is None or old>float(obj@z.x)+1e-6
    if not used:x=z.x
    value=float(obj@x);dual=getattr(z,'mip_dual_bound',None);dual=float(dual) if dual is not None and np.isfinite(dual) else None
    gap=max(0.,-dual-value)/max(1.,abs(value)) if dual is not None else None
    record=dict(label=label,stage=stage+1,objective=key,status=int(z.status),value=value,gap=gap,seconds=time.perf_counter()-now,dual_bound=dual,message=z.message,certificate_kind='HiGHS',used_previous_feasible_incumbent=bool(used))
   logs.append(record)
   mat=vstack([mat,coo_matrix(obj.reshape(1,-1))]).tocsc();lo=np.r_[lo,value-1e-5];hi=np.r_[hi,np.inf]
  np.savez_compressed(folder/(label+'.npz'),x=x,A=self.A,a=self.a,cfg=json.dumps(self.cfg),policy=policy)
  return x,logs

 def evaluate(self,x):
  c=self.cfg;T=self.T;started=np.zeros((2,T));complete=np.zeros((2,T));ontime=np.zeros(2);late=np.zeros(2);wait=0
  for ix,k,r,s,dur in self.jobs:
   n=int(round(x[ix]));started[k,s]+=n;complete[k,s+dur-1]+=n
   if s-r<=c['deadline_h'][k]:ontime[k]+=n
   else:late[k]+=n
   wait+=n*(s-r)
  total=self.A.sum(axis=1);done=started.sum(axis=1);unserved=total-done
  burden=sum(float(self.A[k,r])*(T-r) for k in range(2) for r in range(self.A.shape[1]))+self.scores['burden']@x
  occ=np.array([[sum(x[i] for i in self.occ[k][t])+self.fixed[k,t] for t in range(T)] for k in range(2)])
  q=np.cumsum(np.pad(self.A,((0,0),(0,T-self.A.shape[1])))-started,axis=1)
  assert np.min(q)>-1e-5
  result=dict(arrivals=int(sum(total)),class0_arrivals=int(total[0]),class1_arrivals=int(total[1]),timely_class0=int(round(ontime[0])),timely_class1=int(round(ontime[1])),late_starts=int(round(sum(late))),completed=int(round(sum(done))),unserved=int(round(sum(unserved))),deadline_coverage_pct=100*sum(ontime)/sum(total),weighted_deadline_coverage_pct=100*np.dot(ontime,c['urgency_weights'])/np.dot(total,c['urgency_weights']),queue_patient_hours=float(burden),mean_wait_started_h=wait/max(1,sum(done)),peak_waiting_patients=int(round(q.sum(axis=0).max())),patient_electricity_kwh=float(self.scores['energy']@x),generator_litres=float(self.scores['fuel']@x),spill_kwh=float(sum(x[self.idx['spill_kw']])),generator_kwh=float(sum(x[self.idx['generator_kw']])),all_stages_certified=None)
  hourly=[]
  for t in range(T):
   row=dict(hour=t,grid_available=int(self.a[t]),class0_arrivals=int(self.A[0,t]) if t<self.A.shape[1] else 0,class1_arrivals=int(self.A[1,t]) if t<self.A.shape[1] else 0,class0_starts=int(round(started[0,t])),class1_starts=int(round(started[1,t])),class0_active=float(occ[0,t]),class1_active=float(occ[1,t]),queue=int(round(sum(q[:,t]))),generator_kw=float(sum(x[self.idx['generator_kw'].reshape(-1,T)[:,t]])),generator_on=int(round(sum(x[self.idx['generator_on'].reshape(-1,T)[:,t]]))),fuel_litres=float(x[self.idx['fuel_litres'][t]]),next_fuel_litres=float(x[self.idx['fuel_litres'][t+1]]))
   for name in ['plant','grid','charge_kw','discharge_kw','soc_kwh','spill_kw']:row[name]=float(x[self.idx[name][t]])
   row['next_soc_kwh']=float(x[self.idx['soc_kwh'][t+1]]);hourly.append(row)
  return result,hourly
