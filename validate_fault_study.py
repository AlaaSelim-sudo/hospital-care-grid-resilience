"""Reconstruct service/resource metrics and check source and AC evidence."""
from run_fault_study import *
import hashlib
from collections import defaultdict
def main():
 hashes=json.loads((P/'Source_Hashes.json').read_text())
 for n,d in hashes.items():assert hashlib.sha256((P/Path(n.replace('\\','/'))).read_bytes()).hexdigest()==d,n
 assert hashlib.sha256((P/'figures/county_outage_maps.tikz').read_bytes()).hexdigest()=='8100b59af7d96dc7aee5867d89469aae2cfa795f34e048643d11625772ecd6db'
 snapshot=json.loads((P/'clinical_model/Care_Input_Snapshot.json').read_text())
 nominal=next(r for r in snapshot if float(r['surge'])==1.)
 weights=np.array(nominal['arrivals'])[:,:24];weights=weights/weights.sum()
 evidence=read('Evidence_Daily.csv');cohorts_checked=0
 for case in ['dec2023','apr2024']:
  for h,fips in enumerate(['23001','23011','23019']):
   observed=np.load(P/f'clinical_{case}_H{h+1}.npz')['A'];expected=np.zeros_like(observed)
   er=sorted([r for r in evidence if r['case']==case and r['fips']==fips],key=lambda r:r['date_local'])
   assert len(er)==7
   for day,r in enumerate(er):
    count=int(math.floor(88*float(r['derived_ED_ratio'])+.5));v=(weights*count).ravel();base=np.floor(v).astype(int)
    base[np.argsort(-(v-base),kind='stable')[:count-base.sum()]]+=1
    expected[:,day*24:day*24+24]=base.reshape(2,24)
   assert np.array_equal(observed,expected);cohorts_checked+=1
 checks=[];max_balance=max_soc=max_load=0.;nac=0
 for prefix in ['Fault_','Sensitivity_']:
  results=read(prefix+'Results.csv');starts=read(prefix+'Treatment_Starts.csv');series=read(prefix+'Time_Series.csv');curr=read(prefix+'Current_Loading.csv');nodes=read(prefix+'Node_Voltages.csv')
  for r in results:
   key=r['scenario'];policy=r['policy'];s=scenario(r['case'],repair=int(r['repair_h']),tie_delay=int(r['island_h']),rating=int(r['xfm_kva']),tie_rating=int(r['tie_rating_A']))
   zz=np.load(P/'solutions'/f'{key}_{policy}.npz');occ=np.zeros((3,2,H));occ[:,0,:2]=4;occ[:,1,:1]=8
   admission=np.zeros((3,H));completed=timely=0;served=defaultdict(int)
   for row in starts:
    if row['scenario']!=key or row['policy']!=policy:continue
    h=int(row['hospital'][1:])-1;k=int(row['care_class']);ar=int(row['arrival_hour']);st=int(row['start_hour']);dur=int(row['duration_h']);n=int(row['patients'])
    assert dur==[4,2][k] and ar<=st<=ar+12 and st+dur<=H
    occ[h,k,st:st+dur]+=n
    if k==0:admission[h,st]+=n
    served[h,k,ar]+=n;completed+=n;timely+=n*(st-ar<=[2,8][k])
   assert all(n<=s['A'][h,k,ar] for (h,k,ar),n in served.items())
   assert completed==int(r['completed']) and timely==int(r['timely']) and int(r['arrivals'])==s['A'].sum()
   load=60*zz['plant']+4*occ[:,0]+1.5*occ[:,1]+10*admission
   loaderr=float(np.max(abs(load-zz['load'])));max_load=max(max_load,loaderr);assert loaderr<1e-5
   assert np.max(occ.sum(axis=1)-35*zz['plant'])<1e-5
   assert np.max(2*occ[:,0]+occ[:,1]-20*zz['plant'])<1e-5
   bal=max(float(np.max(abs(zz['p']-(load+zz['charge']-zz['discharge']-zz['generator'])))),float(np.max(abs(zz['q']-(KQ*load-zz['q_support'])))))
   soc=float(np.max(abs(np.diff(zz['energy'],axis=1)-.95*zz['charge']+zz['discharge']/.95)))
   max_balance=max(max_balance,bal);max_soc=max(max_soc,soc);assert bal<1e-5 and soc<1e-5
   assert zz['energy'].min()>=20-1e-5 and zz['energy'].max()<=200+1e-5
   assert np.min(zz['energy'][:,-1])>=160-1e-5 and np.max(abs(zz['energy'][:,0]-200))<1e-5
   assert np.max(zz['charge']*zz['discharge'])<1e-4 and zz['p'].min()>-1e-5
   assert np.max(np.hypot(zz['discharge']-zz['charge'],zz['q_support']))<=110+1e-5
   fuel=np.sum(.30*zz['generator']+1.5*zz['generator_on'],axis=1);assert fuel.max()<=160+1e-5
   assert abs(fuel.sum()-float(r['fuel_litres']))<1e-5
   critical=60*np.sum(1-zz['plant']);assert abs(critical-float(r['critical_unserved_kwh']))<1e-5
   tt=[x for x in series if x['scenario']==key and x['policy']==policy];assert len(tt)==288
   for row in tt:
    assert .5<float(row['min_energized_voltage'])<=float(row['max_energized_voltage'])<1.2
    assert 0<float(row['max_current_ratio'])<3
   uv=sum(int(x['violating_node_phases'])>0 for x in tt)*.25;ol=sum(int(x['overloaded_current_entries'])>0 for x in tt)*.25
   assert abs(uv-float(r['voltage_violation_h']))<1e-8 and abs(ol-float(r['overload_h']))<1e-8
   if policy in ['coordinated','energy_first']:assert uv==ol==0,(key,policy,uv,ol)
   cc=[x for x in curr if x['scenario']==key and x['policy']==policy]
   assert abs(max(float(x['max_current_ratio']) for x in cc)-float(r['max_current_ratio']))<1e-8
   # Replay selected worst cases plus every topology boundary from saved controls.
   f=Network(rating=s['xfm_kva'],tie_rating=s['tie_rating_A'])
   replay={0,47,48,71,72,95,96,143,287,int(float(min(tt,key=lambda x:float(x['min_energized_voltage']))['hour'])*4),int(float(max(tt,key=lambda x:float(x['max_current_ratio']))['hour'])*4)}
   for ix in sorted(replay):
    t=ix//4;stage=s['stages'][t]
    if policy=='island_only' and stage=='backfeed':
     f.topology(stage)
     for h in [1,2]:
      gc.dss.Text.Command(f'edit line.H{h}Supply enabled=no');gc.dss.Text.Command(f'edit load.H{h} enabled=no')
    ac=f.solve_state(stage,background(s['start_hour']+ix*.25),zz['p'][:,t],zz['q'][:,t]);r0=tt[ix]
    assert abs(min(ac['v'][ac['energized']])-float(r0['min_energized_voltage']))<1e-6
    assert abs(max(ac['loading'])-float(r0['max_current_ratio']))<1e-5
    nac+=1
   logs=json.loads((P/'solutions'/f'{key}_{policy}.json').read_text())
   clinical_valid=any(x['stage'] in ['care','energy','clinical_tie_break'] and x['objective'] is not None for x in logs)
   checks.append(dict(scenario=key,policy=policy,timely=timely if clinical_valid else None,completed=completed if clinical_valid else None,clinical_comparison_valid=clinical_valid,critical_deficit_kwh=critical,AC_feasible=uv==ol==0,all_arrivals_completed_energy_ceiling_met=completed==int(s['A'].sum()),solver_time_limits=sum(x['status']==1 for x in logs),max_reported_gap=max([x['gap'] for x in logs if x.get('gap') is not None] or [0])))
 # Verify every plotted cell against raw, quarter-hour values.
 v=defaultdict(list);i=defaultdict(list)
 for row in read('Fault_Node_Voltages.csv'):
  if row['energized']=='True':v[row['case'],row['policy'],row['node'].split('.')[0],int(float(row['hour']))].append(float(row['voltage_pu']))
 for row in read('Fault_Current_Loading.csv'):i[row['case'],row['policy'],row['element'],int(float(row['hour']))].append(float(row['max_current_ratio']))
 cells=read('Plotted_Hourly_Heatmap_Cells.csv')
 for row in cells:
  vals=(v if row['kind']=='minimum_voltage_pu' else i)[row['case'],row['policy'],row['position'],int(row['hour'])]
  if not vals:assert row['value']==''
  else:
   val=min(vals) if row['kind']=='minimum_voltage_pu' else max(vals)
   assert abs(val-float(row['value']))<1e-12
   assert (.90<=val<=1.05) if row['kind']=='minimum_voltage_pu' else (0<=val<=1.20)
 vq=defaultdict(list);iq={}
 for row in read('Fault_Node_Voltages.csv'):
  if row['energized']=='True':vq[row['case'],row['policy'],row['node'].split('.')[0],float(row['hour'])].append(float(row['voltage_pu']))
 for row in read('Fault_Current_Loading.csv'):iq[row['case'],row['policy'],row['element'],float(row['hour'])]=float(row['max_current_ratio'])
 delta_cells=read('Plotted_Hourly_Changes.csv')
 for row in delta_cells:
  h=int(row['hour']);dd=[]
  for k in range(4):
   aa=row['case'],'autonomous',row['position'],h+k*.25;bb=row['case'],'coordinated',row['position'],h+k*.25
   if row['kind']=='paired_voltage_change_pu':
    if vq[aa] and vq[bb]:dd.append(min(vq[bb])-min(vq[aa]))
   elif aa in iq and bb in iq:dd.append(iq[bb]-iq[aa])
  if not dd:assert row['value']==''
  else:
   assert abs(sum(dd)/len(dd)-float(row['value']))<1e-12
   assert abs(float(row['value']))<=float(row['colorbar_abs_limit'])
 report=dict(source_hashes_checked=len(hashes),original_map_unchanged=True,full_week_patient_cohort_arrays_verified=cohorts_checked,schedules_checked=len(checks),selected_AC_replays=nac,max_clinical_load_error_kw=max_load,max_power_balance_error_kw=max_balance,max_SOC_error_kwh=max_soc,heatmap_cells_checked=len(cells),paired_change_cells_checked=len(delta_cells),all_heatmap_scales_shared_no_clipping=True,energization='Source-connected graph; floating disconnected-bus voltages excluded; kV bases recalculated after PCC topology changes',validation_scope='Independent resource/cohort reconstruction and selected AC replay using the same engine and inputs; no field or independent-engine validation',cases=checks)
 (P/'Study_Validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
