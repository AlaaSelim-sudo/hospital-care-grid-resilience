"""Care-constrained grid dispatch with unbalanced OpenDSS verification.

Synthetic benchmark embedding; NOT a hospital feeder or measured voltage study.
One day of the existing assumed arrival profile, plus 24 h clearance/recharge.
"""
from pathlib import Path
import sys,json,csv,math,time,hashlib
P=Path(__file__).resolve().parent;ROOT=P.parents[1]
sys.path[:0]=[str(ROOT/'research/.deps'),str(ROOT/'research/.grid_deps')]
import numpy as np
import scipy
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
import opendssdirect as dss
from care_model_snapshot import EDModel

def write(name,rows):
 if not rows:return
 with (P/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)

class Feeder:
 def __init__(self):
  dss.Text.Command('compile "'+str(P/'feeder/IEEE13Nodeckt.dss')+'"')
  self.original_base=dict(vmin=min(dss.Circuit.AllBusMagPu()),vmax=max(dss.Circuit.AllBusMagPu()))
  dss.Text.Command('set controlmode=off')
  # Explicit scenario modification: lower phases A/C by one tap step so upper
  # voltage is below 1.05. This is a declared utility setting, not a field tap.
  self.taps={'reg1':1.05,'reg2':1.0375,'reg3':1.05}
  for name,tap in self.taps.items():
   dss.Transformers.Name(name);dss.Transformers.Wdg(2);dss.Transformers.Tap(tap)
  self.base=[]
  for name in dss.Loads.AllNames():
   dss.Loads.Name(name);self.base.append((name,dss.Loads.kW(),dss.Loads.kvar()))
  dss.Text.Command('new load.Hospital bus1=675.1.2.3 phases=3 conn=wye kv=4.16 kw=0 kvar=0 model=1 vminpu=.8 vmaxpu=1.2')
  dss.Text.Command('set tolerance=1e-10 maxiterations=100')
  # Declared synthetic one-mile three-phase hospital supply lateral.
  dss.Text.Command('new line.HospitalSupply phases=3 bus1=675.1.2.3 bus2=Hospital.1.2.3 linecode=mtx606 length=1 units=mi')
  dss.Text.Command('edit load.Hospital bus1=Hospital.1.2.3')
  dss.Text.Command('calcvoltagebases')
  self.cache={}
 def solve(self,scale,p,q):
  for name,kw,kvar in self.base:
   dss.Loads.Name(name);dss.Loads.kW(kw*scale);dss.Loads.kvar(kvar*scale)
  dss.Loads.Name('Hospital');dss.Loads.kW(float(p));dss.Loads.kvar(float(q))
  dss.Solution.Solve()
  if not dss.Solution.Converged():raise RuntimeError('OpenDSS did not converge')
  names=dss.Circuit.AllNodeNames();v=np.array(dss.Circuit.AllBusMagPu())
  losses=dss.Circuit.Losses();dss.Circuit.SetActiveBus('Hospital');hv=np.array(dss.Bus.puVmagAngle()[::2])
  return dict(v=v,nodes=names,loss_kw=float(losses[0]/1000),vmin=float(v.min()),vmax=float(v.max()),hospital_vmin=float(hv.min()),weak_node=names[int(v.argmin())],converged=True)
 def linear(self,scale,p0=80.,q0=80.*math.tan(math.acos(.95))):
  key=round(float(scale),8),p0,q0
  if key not in self.cache:
   z=self.solve(scale,p0,q0);vp=self.solve(scale,p0+1,q0)['v'];vq=self.solve(scale,p0,q0+1)['v']
   self.cache[key]=dict(v=z['v'],dp=vp-z['v'],dq=vq-z['v'],p0=p0,q0=q0,nodes=z['nodes'])
  return self.cache[key]

class GridCare(EDModel):
 def __init__(self,arrivals,available,cfg,policy,scales,feeder,pref=None,credit=.0,wait_cap=None):
  super().__init__(arrivals,available,cfg)
  self.policy=policy;T=self.T;old=self.n;self.kq=math.tan(math.acos(.95))
  def var(name,size,low,high,integer=False):
   ids=np.arange(len(self.lb),len(self.lb)+size);self.idx[name]=ids
   self.lb.extend(np.broadcast_to(low,(size,)).tolist());self.ub.extend(np.broadcast_to(high,(size,)).tolist());self.integ.extend([int(integer)]*size);return ids
  qp=var('q_support_pos',T,0,50*np.asarray(available));qn=var('q_support_neg',T,0,50*np.asarray(available))
  reward=var('credited_event_kwh',T,0,300);paid=var('paid_event',T,0,1,True)
  self.n=len(self.lb)
  for key in self.scores:self.scores[key]=np.pad(self.scores[key],(0,self.n-old))
  self.load=[];self.fixed_load=[]
  for t in range(T):
   terms={int(self.idx['plant'][t]):cfg['shared_kw']}
   for k in range(2):
    for ix in self.occ[k][t]:terms[ix]=terms.get(ix,0)+cfg['occupied_bay_kw'][k]
   for ix in self.startA[t]:terms[ix]=terms.get(ix,0)+cfg['diagnostic_kwh_per_class0_start']
   self.load.append(terms);self.fixed_load.append(float(sum(self.fixed[k,t]*cfg['occupied_bay_kw'][k] for k in range(2))))
   # Base case assumes no verified permission for grid-parallel diesel.
   for ix in self.idx['generator_on'].reshape(-1,T)[:,t]:self.ub[ix]=1-int(available[t])
   # Hold an emergency reserve every hour, and restore starting energy by end.
   self.lb[self.idx['soc_kwh'][t]]=max(self.lb[self.idx['soc_kwh'][t]],.5*cfg['battery_kwh'])
   if policy=='care_only':
    if available[t]:self.ub[self.idx['discharge_kw'][t]]=0
    self.ub[qp[t]]=self.ub[qn[t]]=0
   # A conservative inscribed 16-sided inverter capability polygon.
   rating=1.1*cfg['battery_kw']
   for angle in (np.arange(16)+.5)*2*math.pi/16:
    co,si=math.cos(angle),math.sin(angle)
    self.row({int(self.idx['discharge_kw'][t]):co,int(self.idx['charge_kw'][t]):-co,int(qp[t]):si,int(qn[t]):-si},hi=rating*math.cos(math.pi/16))
   eligible=t==16 and all(available[16:19]) and credit>0 and pref is not None
   if not eligible:self.ub[reward[t]]=0;self.ub[paid[t]]=0
   if policy=='coordinated':
    z=feeder.linear(scales[t]);self.add_voltage(t,z,.95005,1.04995)
  # Settle the WHOLE three-hour event. Negative import performance and battery
  # charging reduce net delivery; they cannot be clipped hour by hour.
  if self.ub[reward[16]]>0:
   self.row({int(reward[16]):1,int(paid[16]):-300},hi=0)
   self.row({int(reward[16]):1,int(paid[16]):3000,**{int(self.idx['grid'][t]):1 for t in range(16,19)}},hi=float(sum(pref[16:19]))+3000)
   terms={int(reward[16]):1,int(paid[16]):3000}
   for t in range(16,19):terms[int(self.idx['discharge_kw'][t])]=-1;terms[int(self.idx['charge_kw'][t])]=1
   self.row(terms,hi=3000)
  self.lb[self.idx['soc_kwh'][-1]]=cfg['battery_kwh']
  # Hard care floor: every new visit starts within its class delay target.
  for (k,r),ids in self.bycohort.items():
   timely={ix:1 for ix,kk,rr,s,dur in self.jobs if kk==k and rr==r and s-r<=cfg['deadline_h'][k]}
   self.row(timely,lo=int(self.A[k,r]),hi=int(self.A[k,r]))
  if wait_cap is not None:
   no_start=sum(float(self.A[k,r])*(T-r) for k in range(2) for r in range(self.A.shape[1]))
   self.row({i:float(v) for i,v in enumerate(self.scores['burden']) if v},hi=wait_cap-no_start+1e-5)
  self.cost=np.zeros(self.n);self.prices=np.array([.40 if 16<=t%24<19 else (.12 if t%24<7 or t%24>=22 else .20) for t in range(T)])
  self.cost[self.idx['grid']]=self.prices;self.cost+=1.20*self.scores['fuel']
  self.cost[self.idx['charge_kw']]+=.02;self.cost[self.idx['discharge_kw']]+=.02
  self.cost[qp]+=.0002;self.cost[qn]+=.0002;self.cost[reward]-=credit
  self.credit=credit;self.rebuild()
 def add_voltage(self,t,z,vmin,vmax):
  # q_PCC = availability * (care_load * tan(acos(PF)) - inverter_Q).
  a=float(self.a[t]);g=int(self.idx['grid'][t]);qp=int(self.idx['q_support_pos'][t]);qn=int(self.idx['q_support_neg'][t])
  for v,dp,dq in zip(z['v'],z['dp'],z['dq']):
   coeff={g:float(dp),qp:float(-a*dq),qn:float(a*dq)}
   for ix,value in self.load[t].items():coeff[ix]=coeff.get(ix,0)+float(a*dq*self.kq*value)
   intercept=float(v-dp*z['p0']-dq*z['q0']+a*dq*self.kq*self.fixed_load[t])
   self.row(coeff,lo=vmin-intercept,hi=vmax-intercept)
 def rebuild(self):
  rr=[];cc=[];vv=[]
  for j,terms in enumerate(self.rows):
   for ix,value in terms.items():
    if value:rr.append(j);cc.append(ix);vv.append(value)
  self.M=coo_matrix((vv,(rr,cc)),shape=(len(self.rows),self.n)).tocsc();self.lower=np.array(self.rlo);self.upper=np.array(self.rhi)
 def optimize(self,label,primary_wait=False,seconds=12):
  logs=[]
  if primary_wait:
   z=milp(self.scores['burden'],integrality=self.integ,bounds=Bounds(self.lb,self.ub),constraints=LinearConstraint(self.M,self.lower,self.upper),options={'time_limit':5,'mip_rel_gap':1e-6})
   logs.append(dict(label=label,stage='waiting',status=int(z.status),message=z.message,gap=getattr(z,'mip_gap',None)))
   if z.x is None:return None,logs
   val=float(self.scores['burden']@z.x);self.row({i:float(v) for i,v in enumerate(self.scores['burden']) if v},hi=val+1e-5);self.rebuild()
  z=milp(self.cost,integrality=self.integ,bounds=Bounds(self.lb,self.ub),constraints=LinearConstraint(self.M,self.lower,self.upper),options={'time_limit':seconds,'mip_rel_gap':1e-5})
  logs.append(dict(label=label,stage='cost',status=int(z.status),message=z.message,gap=getattr(z,'mip_gap',None),objective=getattr(z,'fun',None),dual_bound=getattr(z,'mip_dual_bound',None)))
  return z.x,logs

def main():
 feeder=Feeder();T=48
 shape=np.array([.96,.95,.94,.94,.95,.96,.97,.98,.99,1.,1.,1.,1.,1.,1.,1.01,1.02,1.02,1.02,1.01,1.,.99,.98,.97])
 scales=np.tile(shape-.0250,2)
 original=json.loads((P/'Care_Input_Snapshot.json').read_text())
 configs=[]
 for mode in ['connected','county_outage_replay']:
  for inp in original:
   cfg=dict(inp['configuration']);A=np.array(inp['arrivals'])[:,:24]
   available=np.ones(T,int)
   if mode=='county_outage_replay':available[:24]=np.array(inp['available'])[:24]
   configs.append((mode,inp,A,available,cfg))
 results=[];hourly=[];node_rows=[];logs=[];input_rows=[]
 for mode,inp,A,a,cfg in configs:
  surge=float(inp['surge']);scenario=f'{mode}_s{surge}'
  ref_p=None;wait_cap=None
  for policy in ['care_only','tariff_only','coordinated']:
   label=scenario+'_'+policy
   # Zero assumed local-contract revenue for actual December/April-type storms.
   # Connected case is a distinct proposed voltage-service contract, not NH DR.
   credit=.25 if mode=='connected' and policy!='care_only' else 0.
   model=GridCare(A,a,cfg,policy,scales,feeder,ref_p,credit,wait_cap)
   print('SOLVE',label,model.n,len(model.rows),flush=True)
   x,ll=model.optimize(label,primary_wait=policy=='care_only');logs+=ll
   if x is None:
    results.append(dict(scenario=scenario,policy=policy,status='NO_FEASIBLE_INCUMBENT',solver_status=ll[-1]['status']));continue
   for iteration in range(3):
    _,hs=model.evaluate(x);replay=[]
    for t,h in enumerate(hs):
     load=model.fixed_load[t]+sum(v*x[ix] for ix,v in model.load[t].items())
     qinv=float(x[model.idx['q_support_pos'][t]]-x[model.idx['q_support_neg'][t]])
     qnet=float(a[t])*(model.kq*load-qinv)
     zz=feeder.solve(scales[t],h['grid'],qnet);replay.append((zz,load,qnet,qinv))
    bad=[t for t,(z,*_) in enumerate(replay) if z['vmin']<.95-1e-6 or z['vmax']>1.05+1e-6]
    if policy!='coordinated' or not bad:break
    for t in bad:
     z,load,qnet,qinv=replay[t];lin=feeder.linear(scales[t],float(hs[t]['grid']),qnet);model.add_voltage(t,lin,.95005,1.04995)
    model.rebuild();x,ll=model.optimize(label+'_refine'+str(iteration+1),seconds=15);logs+=ll
    if x is None:break
   if x is None:
    results.append(dict(scenario=scenario,policy=policy,status='REFINEMENT_NO_INCUMBENT'));continue
   # Recompute after the final refinement so no stale AC replay can be saved.
   ev,hs=model.evaluate(x);replay=[]
   for t,h in enumerate(hs):
    load=model.fixed_load[t]+sum(v*x[ix] for ix,v in model.load[t].items())
    qinv=float(x[model.idx['q_support_pos'][t]]-x[model.idx['q_support_neg'][t]])
    qnet=float(a[t])*(model.kq*load-qinv)
    replay.append((feeder.solve(scales[t],h['grid'],qnet),load,qnet,qinv))
   if policy=='coordinated' and any(z['vmin']<.95 or z['vmax']>1.05 for z,*_ in replay):
    results.append(dict(scenario=scenario,policy=policy,status='AC_VERIFICATION_FAILED'));continue
   violation=max(float(np.max(model.M@x-model.upper)),float(np.max(model.lower-model.M@x)),float(np.max(np.array(model.lb)-x)),float(np.max(x-np.array(model.ub))),0.)
   assert violation<2e-5
   assert np.max(np.abs(x[np.array(model.integ,dtype=bool)]-np.round(x[np.array(model.integ,dtype=bool)])))<2e-5
   assert ev['unserved']==0 and ev['late_starts']==0
   if policy=='care_only':ref_p=np.array([h['grid'] for h in hs]);wait_cap=ev['queue_patient_hours']
   assert wait_cap is None or ev['queue_patient_hours']<=wait_cap+1e-4
   e_cost=sum(model.prices[t]*h['grid'] for t,h in enumerate(hs));fuel_cost=1.2*ev['generator_litres'];wear=.02*sum(h['charge_kw']+h['discharge_kw'] for h in hs)
   credits=credit*float(sum(x[model.idx['credited_event_kwh']]))
   uv=sum(z['vmin']<.95-1e-6 for z,*_ in replay);uv_area=sum(float(np.maximum(0,.95-z['v']).sum()) for z,*_ in replay)
   row=dict(scenario=scenario,mode=mode,surge=surge,policy=policy,status='FEASIBLE_SAVED_SCHEDULE',arrivals=ev['arrivals'],timely_pct=ev['deadline_coverage_pct'],waiting_patient_h=ev['queue_patient_hours'],grid_kwh=sum(h['grid'] for h in hs),peak_grid_kw=max(h['grid'] for h in hs),diesel_litres=ev['generator_litres'],grid_cost_dollars=e_cost,fuel_cost_dollars=fuel_cost,battery_throughput_cost_dollars=wear,hypothetical_local_credit_dollars=credits,net_variable_cost_dollars=e_cost+fuel_cost+wear-credits,min_feeder_voltage_pu=min(z['vmin'] for z,*_ in replay),min_hospital_voltage_pu=min(z['hospital_vmin'] for z,*_ in replay),max_feeder_voltage_pu=max(z['vmax'] for z,*_ in replay),undervoltage_hours=uv,undervoltage_node_pu_hours=uv_area,losses_kwh=sum(z['loss_kw'] for z,*_ in replay),max_linear_constraint_residual=violation,solver_status=ll[-1]['status'],cost_gap=ll[-1].get('gap'),AC_refinements=iteration,final_soc_kwh=hs[-1]['next_soc_kwh'])
   results.append(row)
   for t,(h,(z,load,qnet,qinv)) in enumerate(zip(hs,replay)):
    hourly.append(dict(scenario=scenario,policy=policy,**h,care_load_kw=load,net_q_kvar=qnet,inverter_q_kvar=qinv,background_multiplier=float(scales[t]),min_feeder_pu=z['vmin'],max_feeder_pu=z['vmax'],hospital_min_pu=z['hospital_vmin'],loss_kw=z['loss_kw'],credited_event_kwh=float(x[model.idx['credited_event_kwh'][t]])))
    for name,v in zip(z['nodes'],z['v']):node_rows.append(dict(scenario=scenario,policy=policy,hour=t,node=name,voltage_pu=float(v)))
   np.savez_compressed(P/(label+'.npz'),x=x,A=A,a=a,cfg=json.dumps(cfg))
   input_rows.append(dict(label=label,configuration=cfg,arrivals=A.tolist(),available=a.tolist(),background_multiplier=scales.tolist(),credit_dollars_per_kwh=credit,reference_grid_kw=ref_p.tolist() if ref_p is not None else None,waiting_cap_patient_h=wait_cap))
   write('Pilot_Results.csv',results);write('Pilot_Hourly.csv',hourly);write('Node_Voltages.csv',node_rows)
   (P/'Solver_Log.json').write_text(json.dumps(logs,indent=2));(P/'Pilot_Inputs.json').write_text(json.dumps(input_rows,indent=2))
   print('RESULT',row,flush=True)
 # Report a sweep that includes states beyond the chosen manageable case.
 sweep=[]
 for scale in [1.,1.01,1.02,1.03,1.04]:
  for kw in [0,60,90,120,150]:
   for qsupport in [0,50]:
    z=feeder.solve(scale,kw,kw*model.kq-qsupport)
    sweep.append(dict(background_multiplier=scale,hospital_import_kw=kw,q_support_kvar=qsupport,min_voltage_pu=z['vmin'],max_voltage_pu=z['vmax'],weak_node=z['weak_node']))
 write('Voltage_Stress_Sweep.csv',sweep)
 metadata=dict(engine=dss.Basic.Version(),scipy=scipy.__version__,feeder='Modified IEEE 13-node benchmark; hospital balanced PQ attachment via a one-mile mtx606 lateral from bus 675; common background offset -0.0250',original_benchmark=feeder.original_base,scenario_taps=feeder.taps,voltage_band=[.95,1.05],linear_guard_band=[.95005,1.04995],hospital_pf=.95,inverter_kva=110,q_limit_kvar=50,arrival_hours=24,clearance_hours=24,lateral_miles=1,background_offset=-.0250,source_care_profile='First 24 h of existing assumed April profiles; not observed storm ED arrivals',costs=dict(grid_night_dollars_per_kwh=.12,grid_day=.20,grid_peak=.40,diesel_dollars_per_litre=1.20,battery_dollars_per_kwh_throughput=.02,proposed_connected_local_credit_dollars_per_kwh=.25),incentive_scope='Connected-case proposed local contract only; not existing NH/VT utility payment. Winter outage replay credits zero.',parallel_diesel=False,reserve='SOC >=50% every hour; terminal SOC equals initial 100%',not_modeled='Measured feeder, campus HVAC, heat storage, real program eligibility, current tariffs, inverter losses, thermal constraints, dynamic ride-through, islanding transients, reconnection, generator reactive capacity',cases=len(results))
 (P/'Pilot_Metadata.json').write_text(json.dumps(metadata,indent=2));print('DONE',len(results),'cases',flush=True)
if __name__=='__main__':main()
