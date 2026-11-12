"""Make solver and clinical interpretation explicit in the downloadable tables."""
from run_fault_study import *
def main():
 for prefix in ['Fault_','Sensitivity_']:
  rows=read(prefix+'Results.csv')
  for r in rows:
   log=json.loads((P/'solutions'/f"{r['scenario']}_{r['policy']}.json").read_text())
   patient=[s for s in log if s['stage'] in ['care','energy','clinical_tie_break']]
   valid=any(s['objective'] is not None for s in patient)
   r['clinical_comparison_valid']=valid
   r['patient_result_scope']='Accepted clinical incumbent; remaining gaps in solver logs' if valid else 'Critical-service feasibility solution only; patient optimization unresolved; not a clinical comparison result'
   r['critical_service_stage_gap']=log[0].get('gap')
   r['patient_stage_gaps']=json.dumps({s['stage']:s.get('gap') for s in patient})
   r['all_arrivals_completed_energy_ceiling_met']=int(r['completed'])==int(r['arrivals'])
  write(prefix+'Results.csv',rows)
 evidence=read('Evidence_Daily.csv');inputs=[]
 for case in ['dec2023','apr2024']:
  s=scenario(case)
  for h,(fips,county) in enumerate([('23001','Androscoggin'),('23011','Kennebec'),('23019','Penobscot')]):
   for t in range(48):
    local=datetime.fromisoformat(s['window_start'])+timedelta(hours=t)
    row=next(r for r in evidence if r['case']==case and r['fips']==fips and r['date_local']==str(local.date()))
    inputs.append(dict(case=case,hospital=f'H{h+1}',county=county,relative_hour=t,local_time=local.isoformat(),local_zone='EST' if case=='dec2023' else 'EDT',observed_county_ED=row['observed_ED'],estimated_county_weekday_ED=row['estimated_weekday_ED'],derived_county_demand_ratio=row['derived_ED_ratio'],modeled_class0_arrivals=int(s['A'][h,0,t]),modeled_class1_arrivals=int(s['A'][h,1,t]),scope='County totals observed; hospital transfer, hourly pattern and acuity split assumed'))
 write('Fault_Demand_Inputs.csv',inputs)
 print('Annotated all rows; critical-only candidate counts are explicitly excluded from clinical comparisons.')
if __name__=='__main__':main()
