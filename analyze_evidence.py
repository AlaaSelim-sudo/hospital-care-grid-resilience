from pathlib import Path
import csv,json,math,statistics,hashlib
from datetime import datetime,date,timedelta,timezone
from collections import defaultdict
P=Path(__file__).resolve().parent
def read(n):return list(csv.DictReader((P/n).open(encoding='utf-8-sig')))
def write(n,rr):
 with (P/n).open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def corr(x,y):return statistics.correlation(x,y) if len(x)>=5 and len(set(x))>1 and len(set(y))>1 else None
def main():
 observed=read('Observed_County_Daily_ED.csv');ed={(r['case'],r['fips'],r['date_local']):int(r['observed_ED_visits']) for r in observed}
 rawed={}
 for src in json.loads((P/'ED_Download_Provenance.json').read_text()):
  path=P/src['file'].replace('\\','/');bb=path.read_bytes();assert hashlib.sha256(bb).hexdigest()==src['sha256'];rawed[src['file']]=list(csv.reader(bb.decode('utf-16').splitlines(),delimiter='\t'))
 for r in observed:
  row=rawed[r['source_file']][int(r['ED_source_row'])-1];assert row[10]=='Visit Count'
  values=[float(v.replace(',','')) for v in row[11:] if v!=''];assert len(values)==1 and values[0]==int(r['observed_ED_visits'])
  assert str(datetime.strptime(row[4],'%b %d, %Y').date())==r['date_local']
 events={'dec2023':('2023-12-18',-5),'apr2024':('2024-04-03',-4)};daily=[];correlations=[];summary=[];checks=[]
 for case,(startstr,tzoff) in events.items():
  start=date.fromisoformat(startstr);raw=read('data/'+case+'_EAGLEI_NH_VT_ME.csv');prov=json.loads((P/'data'/f'{case}_extraction_provenance.json').read_text())
  assert hashlib.sha256((P/'data'/f'{case}_EAGLEI_NH_VT_ME.csv').read_bytes()).hexdigest()==prov['regional_extract_sha256']
  byday=defaultdict(list)
  for r in raw:
   if not r['fips_code'].startswith('23'):continue
   t=datetime.fromisoformat(r['run_start_time']).replace(tzinfo=timezone.utc).astimezone(timezone(timedelta(hours=tzoff)))
   byday[r['fips_code'],t.date()].append(float(r.get('customers_out',r.get('sum'))))
  for fips in sorted({r['fips'] for r in observed if r['case']==case}):
   county=next(r['county'] for r in observed if r['fips']==fips)
   pre=[start-timedelta(days=j) for j in range(21,0,-1)]
   baseline={d:statistics.mean(ed[case,fips,str(dt)] for dt in pre if dt.weekday()==d) for d in range(7)}
   baseline14={d:statistics.mean(ed[case,fips,str(dt)] for dt in pre[-14:] if dt.weekday()==d) for d in range(7)}
   pairs=[]
   for day in range(7):
    dt=start+timedelta(days=day);out=byday[fips,dt];base=baseline[dt.weekday()];ob=ed[case,fips,str(dt)];H=.25*sum(out) if out else None
    nextdt=dt+timedelta(days=1);nextob=ed[case,fips,str(nextdt)]
    row=dict(case=case,fips=fips,county=county,date_local=str(dt),observed_ED=ob,estimated_weekday_ED=base,derived_ED_ratio=ob/base,derived_ED_change_pct=100*(ob/base-1),recorded_customer_outage_h='' if H is None else H,snapshots=len(out),coverage=len(out)/96,derived_next_day_ED_change_pct=100*(nextob/baseline[nextdt.weekday()]-1),derived_ED_change_14day_pct=100*(ob/baseline14[dt.weekday()]-1),derived_next_day_ED_change_14day_pct=100*(nextob/baseline14[nextdt.weekday()]-1))
    daily.append(row)
    if H is not None:pairs.append(row)
   for lag,ykey,y14 in [(0,'derived_ED_change_pct','derived_ED_change_14day_pct'),(1,'derived_next_day_ED_change_pct','derived_next_day_ED_change_14day_pct')]:
    for threshold in [0,.75,.9]:
     subset=[r for r in pairs if r['coverage']>=threshold];x=[r['recorded_customer_outage_h'] for r in subset];y=[r[ykey] for r in subset]
     rv=corr(x,y);r14=corr(x,[r[y14] for r in subset])
     correlations.append(dict(case=case,fips=fips,county=county,lag_days=lag,min_snapshot_coverage=threshold,n_days=len(subset),pearson_r='' if rv is None else rv,r_14day_baseline='' if r14 is None else r14,method='Recorded customer-hours vs estimated weekday-adjusted ED change; exploratory, no causal or significance claim'))
  for lag in [0,1]:
   for th in [0,.75,.9]:
    vals=[float(r['pearson_r']) for r in correlations if r['case']==case and r['lag_days']==lag and r['min_snapshot_coverage']==th and r['pearson_r']!='']
    summary.append(dict(case=case,lag_days=lag,min_snapshot_coverage=th,estimable_counties=len(vals),positive_counties=sum(v>0 for v in vals),median_r=statistics.median(vals) if vals else '',min_r=min(vals) if vals else '',max_r=max(vals) if vals else ''))
  checks.append(dict(case=case,raw_rows=len(raw),hash_verified=True))
 # Verify previously published county primary coefficients from source recomputation.
 old=read('County_Event_Results.csv')
 for r in correlations:
  if r['min_snapshot_coverage']!=0:continue
  o=next(x for x in old if x['case']==r['case'] and x['fips']==r['fips'])
  key='pearson_same_day_baseline_adjusted_ED_r' if r['lag_days']==0 else 'pearson_lag1_baseline_adjusted_ED_r'
  if r['pearson_r']!='':assert abs(float(r['pearson_r'])-float(o[key]))<1e-10
 write('Evidence_Daily.csv',daily);write('Evidence_Correlations.csv',correlations);write('Evidence_Summary.csv',summary)
 (P/'Evidence_Validation.json').write_text(json.dumps(dict(county_event_days=len(daily),observed_ED_source_rows=len(observed),all_ED_source_counts_rechecked_against_original_CSV_records=True,ED_source_file_hashes_verified=len(rawed),recomputed_primary_coefficients_match_prior=True,timezone_alignment='UTC to EST/EDT local calendar date',missing_outages_imputed=False,source_extracts=checks),indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
