"""Separate electrical fault snapshots; no protection/transient validation claimed."""
from run_fault_study import *
from fault_network import dss,LINE_RATINGS
def main():
 rows=[]
 for case in ['dec2023','apr2024']:
  s=scenario(case);t=s['fault_hour'];p=s['ref'][:,t].copy();q=p*KQ
  for typ,ph,bus,r in [('SLG',1,'F.3',.5),('3PG',3,'F.1.2.3',.05)]:
   f=Network(rating=s['xfm_kva']);pre=f.solve_state('normal',background(s['start_hour']+t),p,q)
   # Split the damaged section: the fault is inside the isolated line, not on
   # a bus that would be re-energized by the restoration tie.
   dss.Text.Command('edit line.670671 bus2=F.1.2.3 length=666.5 units=ft')
   dss.Text.Command('new line.FaultTail phases=3 bus1=F.1.2.3 bus2=671.1.2.3 linecode=mtx601 length=666.5 units=ft')
   dss.Text.Command('calcvoltagebases')
   ground='F.0' if ph==1 else 'F.0.0.0'
   dss.Text.Command(f'new fault.StormFault phases={ph} bus1={bus} bus2={ground} r={r}')
   dss.Text.Command('set loadmodel=admittance');dss.Solution.SolveDirect()
   assert dss.Solution.Converged()
   dss.Circuit.SetActiveElement('Fault.StormFault');fi=max(dss.CktElement.CurrentsMagAng()[::2]);v=np.array(dss.Circuit.AllBusMagPu());names=dss.Circuit.AllNodeNames()
   hh={f'H{h+1}_minimum_pu':float(min(v[j] for j,n in enumerate(names) if n.split('.')[0]==f'h{h+1}')) for h in range(3)}
   rows.append(dict(case=case,fault_type=typ,fault_bus='F: midpoint of 670--671',fault_resistance_ohm=r,fault_current_A=fi,recorded_county_peak_UTC=s['fault_UTC'],assumed_trip_s=.1,**hh))
 write('Fault_Snapshot_Results.csv',rows)
 config=dict(observed=['Maine residence-county daily ED counts','EAGLE-I county outage timestamps and counts','original NH/VT/ME exposure/geometry'],assumed=['fault at midpoint F of line 670--671; isolate both ends','0.5-ohm phase-C ground and 0.05-ohm 3-phase-ground fault snapshots','0.1-second ideal protective isolation; not a relay simulation','6-h isolation, then temporary tie; physical repair at 12 h','normally-open 633--680 tie, 0.05 mile, mtx601, assumed 500 A','H3 connected behind 634, .48 kV, .005-mile lateral','XFM1 600 kVA capacity scenario; original 500-kVA audit retained separately','634 phase load balanced at unchanged total 400 kW / 290 kvar, common to all policies','fixed regulator taps 1.0495/1.0375/1.0495','three synthetic 88-visits/day hospitals, county daily ratios; no hospital-specific observations','200 kWh / 100 kW / 110 kVA batteries, 20 kWh floor, 200 kWh initial, 160 kWh terminal','100 kW generator per site, 10-kW minimum, 160 litres initial fuel, no refill','0.30 L/kWh + 1.5 L/h when on, assumed $1.20/L','onsite island-forming capability, ideal transfer and parallel operation permission; no dynamics/protection validation'],assumed_line_current_limits_A=LINE_RATINGS,scope='Quasi-static fault snapshots plus multi-hour post-trip care/energy scheduling and 15-minute AC replay. A sustained fault is never dispatched through. County outages do not identify actual damaged feeders or hospital interruption durations.',protocol=['share critical load, SOC, fuel and flexibility bounds','utility announces damaged line, assumed isolation/repair times and common tie opportunity','co-optimize nonpreemptive care starts, generator fuel, battery P/Q and grid import under V/I limits','check every energized phase and winding with OpenDSS','hold unsafe connection or tighten/re-solve; do not count unsafe schedules as guaranteed care'])
 (P/'Study_Assumptions.json').write_text(json.dumps(config,indent=2));print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
