"""Simultaneous hospital surge: fixed-care, multi-PCC dispatch and AC replay.
All hospital locations, equipment and profiles are benchmark assumptions.
"""
from pathlib import Path
import sys,json,csv,math,shutil,time
P=Path(__file__).resolve().parent;ROOT=P.parents[1];BASE=P/'clinical_model'
sys.path[:0]=[str(ROOT/'research/.deps'),str(ROOT/'research/.grid_deps'),str(P)]
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
import opendssdirect as dss
DT=.25;T=192;NH=3;KQ=math.tan(math.acos(.95))
def write(name,rows):
 if not rows:return
 with (P/name).open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
class Feeder:
 def __init__(self):
  dss.Text.Command('compile "'+str(P/'feeder/IEEE13Nodeckt.dss')+'"');dss.Text.Command('set controlmode=off')
  for n,tap in zip(['reg1','reg2','reg3'],[1.05,1.0375,1.05]):
   dss.Transformers.Name(n);dss.Transformers.Wdg(2);dss.Transformers.Tap(tap)
  self.base=[]
  for n in dss.Loads.AllNames():
   dss.Loads.Name(n);self.base.append((n,dss.Loads.kW(),dss.Loads.kvar()))
  for i,(bus,length) in enumerate([('675',1.),('680',.75),('633',.5)]):
   dss.Text.Command(f'new line.H{i+1}Supply phases=3 bus1={bus}.1.2.3 bus2=H{i+1}.1.2.3 linecode=mtx606 length={length} units=mi')
   dss.Text.Command(f'new load.H{i+1} bus1=H{i+1}.1.2.3 phases=3 conn=wye kv=4.16 kw=0 kvar=0 model=1 vminpu=.8 vmaxpu=1.2')
  dss.Text.Command('calcvoltagebases');dss.Text.Command('set tolerance=1e-10 maxiterations=100')
  self.cache={}
 def solve(self,scale,p,q):
  for n,kw,kvar in self.base:
   dss.Loads.Name(n);dss.Loads.kW(kw*scale);dss.Loads.kvar(kvar*scale)
  for i in range(NH):
   dss.Loads.Name(f'H{i+1}');dss.Loads.kW(float(p[i]));dss.Loads.kvar(float(q[i]))
  dss.Solution.Solve();assert dss.Solution.Converged()
  v=np.asarray(dss.Circuit.AllBusMagPu());names=dss.Circuit.AllNodeNames()
  hp=[]
  for i in range(NH):
   active=dss.Circuit.SetActiveBus(f'H{i+1}')
   hp.append(min(dss.Bus.puVmagAngle()[::2]) if active>=0 else 0.)
  return dict(v=v,names=names,hospital_v=hp,loss_kw=dss.Circuit.Losses()[0]/1000)
 def linear(self,scale,p,q):
  key=(round(float(scale),9),tuple(np.round(p,7)),tuple(np.round(q,7)))
  if key not in self.cache:
   z=self.solve(scale,p,q);dp=[];dq=[]
   for h in range(NH):
    pp=np.asarray(p).copy();pp[h]+=1;dp.append(self.solve(scale,pp,q)['v']-z['v'])
    qq=np.asarray(q).copy();qq[h]+=1;dq.append(self.solve(scale,p,qq)['v']-z['v'])
   self.cache[key]=(z['v'],np.asarray(dp).T,np.asarray(dq).T)
  return self.cache[key]
def dispatch(load,scales,f,policy,energy=200.,power=100.,controls='PQ',reserve_fractions=None):
 # Fixed care decisions: no load shedding, transfers or treatment rescheduling.
 idx={};lb=[];ub=[];integ=[]
 def var(n,shape,low,high,integer=0):
  start=len(lb);size=int(np.prod(shape));ix=np.arange(start,start+size).reshape(shape);idx[n]=ix
  lb.extend(np.broadcast_to(low,shape).ravel());ub.extend(np.broadcast_to(high,shape).ravel());integ.extend([integer]*size);return ix
 c=var('charge',(NH,T),0,power);d=var('discharge',(NH,T),0,power);qp=var('qpos',(NH,T),0,.5*power);qn=var('qneg',(NH,T),0,.5*power)
 e=var('energy',(NH,T+1),energy*.5 if reserve_fractions is None else energy*np.asarray(reserve_fractions),energy);binary=var('mode',(NH,T),0,1,1)
 for h in range(NH):lb[e[h,0]]=ub[e[h,0]]=energy;lb[e[h,-1]]=ub[e[h,-1]]=energy
 slack=var('voltage_slack',(T,),0,.2 if policy in ['coordinated','local_voltage'] else 0)
 n=len(lb);rows=[];rlo=[];rhi=[]
 def row(v,lo=-np.inf,hi=np.inf):rows.append(v);rlo.append(lo);rhi.append(hi)
 for h in range(NH):
  for t in range(T):
   row({int(e[h,t+1]):1,int(e[h,t]):-1,int(c[h,t]):-DT*.95,int(d[h,t]):DT/.95},0,0)
   row({int(c[h,t]):1,int(binary[h,t]):-power},hi=0);row({int(d[h,t]):1,int(binary[h,t]):power},hi=power)
   row({int(d[h,t]):1,int(c[h,t]):-1},hi=load[h,t]) # no export
   for angle in (np.arange(16)+.5)*2*np.pi/16:
    co,si=np.cos(angle),np.sin(angle)
    row({int(d[h,t]):co,int(c[h,t]):-co,int(qp[h,t]):si,int(qn[h,t]):-si},hi=1.1*power*np.cos(np.pi/16))
 if policy=='independent':
  for ix in np.r_[qp.ravel(),qn.ravel()]:ub[ix]=0
 if controls=='P':
  for ix in np.r_[qp.ravel(),qn.ravel()]:ub[ix]=0
 if controls=='Q':
  for ix in np.r_[c.ravel(),d.ravel()]:ub[ix]=0
 if policy in ['coordinated','local_voltage']:
  for t in range(T):
   v,dp,dq=f.linear(scales[t],load[:,t],KQ*load[:,t])
   for j in range(len(v)):
    node=dss.Circuit.AllNodeNames()[j]
    if policy=='local_voltage' and node.split('.')[0] not in ['h1','h2','h3']:continue
    terms={int(slack[t]):1};terms_hi={int(slack[t]):-1}
    for h in range(NH):
     if policy=='local_voltage' and node.split('.')[0]!=f'h{h+1}':continue
     for ix,val in [(c[h,t],dp[j,h]),(d[h,t],-dp[j,h]),(qp[h,t],-dq[j,h]),(qn[h,t],dq[j,h])]:terms[int(ix)]=float(val);terms_hi[int(ix)]=float(val)
    row(terms,lo=.9502-v[j]);row(terms_hi,hi=1.04995-v[j])
 def matrix():
  ir=[];ic=[];vs=[]
  for i,r in enumerate(rows):
   for j,v in r.items():ir.append(i);ic.append(j);vs.append(v)
  return coo_matrix((vs,(ir,ic)),shape=(len(rows),n)).tocsc()
 M=matrix();logs=[]
 if policy in ['coordinated','local_voltage']:
  objective=np.zeros(n);objective[slack]=DT
  hardub=np.asarray(ub).copy();hardub[slack]=0
  z=milp(np.zeros(n),integrality=integ,bounds=Bounds(lb,hardub),constraints=LinearConstraint(M,rlo,rhi),options={'time_limit':20,'mip_rel_gap':1e-5})
  logs.append(dict(stage='hard_voltage_feasibility',status=int(z.status),gap=getattr(z,'mip_gap',None),objective=z.fun))
  if z.x is not None:
   for ix in slack:ub[ix]=0
  else:
   z=milp(objective,integrality=integ,bounds=Bounds(lb,ub),constraints=LinearConstraint(M,rlo,rhi),options={'time_limit':40,'mip_rel_gap':1e-5})
   logs.append(dict(stage='minimum_guard_violation',status=int(z.status),gap=getattr(z,'mip_gap',None),objective=z.fun));assert z.x is not None
   row({int(ix):DT for ix in slack},hi=float(z.fun)+1e-8);M=matrix()
 prices=np.array([.4 if 16<=int(t*DT)%24<19 else (.12 if int(t*DT)%24<7 or int(t*DT)%24>=22 else .2) for t in range(T)])
 cost=np.zeros(n)
 for h in range(NH):cost[c[h]]=DT*(prices+.02);cost[d[h]]=DT*(-prices+.02);cost[qp[h]]=cost[qn[h]]=DT*.0002
 z=milp(cost,integrality=integ,bounds=Bounds(lb,ub),constraints=LinearConstraint(M,rlo,rhi),options={'time_limit':45,'mip_rel_gap':1e-5})
 assert z.x is not None;logs.append(dict(stage='cost',status=int(z.status),gap=getattr(z,'mip_gap',None),objective=z.fun))
 x=z.x;res=max(np.max(M@x-np.asarray(rhi)),np.max(np.asarray(rlo)-M@x),np.max(np.asarray(lb)-x),np.max(x-np.asarray(ub)),0.)
 assert res<2e-5
 pp=load+x[c]-x[d];qq=KQ*load-x[qp]+x[qn]
 return dict(p=pp,q=qq,charge=x[c],discharge=x[d],q_support=x[qp]-x[qn],energy=x[e],cost=float(np.sum(pp*prices)*DT+.02*np.sum(x[c]+x[d])*DT),slack=x[slack],residual=float(res),logs=logs,x=x,idx=idx)
