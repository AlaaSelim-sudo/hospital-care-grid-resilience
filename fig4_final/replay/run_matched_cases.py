import argparse,json,csv,sys
from pathlib import Path
root=Path.cwd()
sys.path.insert(0,str(root))
from severe_stress import settings_for,add_island_cut
from compound_stress import solve
from run_fault_study import Model,evaluate,linear_data,np
parser=argparse.ArgumentParser()
parser.add_argument('event')
parser.add_argument('case')
parser.add_argument('policy')
args=parser.parse_args()
out=root/'matched_case_runs'
out.mkdir(exist_ok=True)
stem=out/f'{args.event}_{args.case}_{args.policy}'
s,added=settings_for(args.event,args.case)
assert (s['tie_delay_h'],s['repair_h'],s['fuel_limit_l'],s['backfeed_background_fraction'])==(18,27,(320,320,160),.95)
other,_=settings_for(args.event,'long_no_surge' if args.case=='long_reserve' else 'long_reserve')
assert all(np.array_equal(s[k],other[k]) for k in s if k!='A')
joint=args.policy=='coordinated'
m=Model(s,args.policy,linear_data(s),tighten_v=.0003 if joint else 0,tighten_i=.012 if joint else 0)
add_island_cut(m,18,320)
x,logs=solve(m,45,stem.with_suffix('.log.json'),seed=None,stages=('essential','urgent'),prime=True)
assert x is not None
metrics,trace,*_=evaluate(s,m,x)
assert metrics['linear_residual']<1e-4
np.savez_compressed(stem.with_suffix('.npz'),x=x)
with stem.with_suffix('.replay.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(trace[0]))
    w.writeheader()
    w.writerows(trace)
metrics.update(case=args.case,event=args.event,policy=args.policy,extra_arrivals=added,stages=logs)
stem.with_suffix('.result.json').write_text(json.dumps(metrics,indent=2))
print(args,metrics['min_voltage'],metrics['voltage_violation_h'],[(r['stage'],r['certified'],r['gap']) for r in logs],flush=True)
