# Hospital care and grid resilience under feeder faults

This study schedules emergency care and backup power together during a feeder outage and staged restoration. It follows two storm scenarios, December 2023 and April 2024, on a modified IEEE 13-node feeder with three modeled hospitals.

## 1. Set up and run the study

Download [the study files and datasets](reproduction.zip) and extract them. Open PowerShell in the extracted directory, where `reproduce.py` and `requirements.txt` are located. Use Python 3.12.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe reproduce.py --skip-pdf
```

On macOS or Linux, use `python3.12 -m venv .venv` and replace `.\.venv\Scripts\python.exe` with `./.venv/bin/python`. The full run executes the stages below and writes CSVs and TikZ files into the extracted directory.

## 2. Check the observations and build the scenarios

Run `python analyze_evidence.py` using the environment above. It reads archived EAGLE-I outage snapshots for NH, VT, and ME and Maine's public county ED counts. It aligns records by local calendar date and compares ED visits with a weekday baseline from the preceding three weeks. Inspect `Observed_County_Daily_ED.csv` for the source counts and `Evidence_Daily.csv` for the derived ratios.

In [`run_fault_study.py`](run_fault_study.py), start with `scenario()`. It loads the modeled cohorts and places the simulated fault at the recorded Kennebec County outage peak. The model uses county ED ratios to scale an assumed hospital arrival profile. The hospital equipment, hourly arrivals, fault location, and switching times are scenario inputs. In [`severe_stress.py`](severe_stress.py), `settings_for()` selects the longer outages, surge, fuel, and staged-pickup settings.

## 3. Follow the optimization model

1. Read `Model` in [`run_fault_study.py`](run_fault_study.py). It assigns treatment starts to arrival cohorts, tracks patients still in treatment, and converts treatment and essential services into electrical demand.
2. Follow the same model into the generator, battery, fuel, and grid-import constraints. `linear_data()` uses OpenDSS to estimate how changes in hospital imports affect phase voltage and equipment loading. [`fault_network.py`](fault_network.py) implements the feeder and switching states.
3. Read `run()` in [`severe_stress.py`](severe_stress.py) to see the objective order. For the matched 18-hour comparison, both policies first prioritize essential service. Care-first then maximizes on-time urgent starts. Energy-first maximizes treatment energy and then maximizes on-time urgent starts at the fixed energy outcome.

## 4. Run and check the cases

Use the virtual environment's Python executable for these commands.

```sh
python run_fault_study.py
python compound_stress.py --replay
python severe_stress.py --replay
python severe_stress.py --reserve --replay
python severe_stress.py --certificate --replay
```

These commands check the stored decisions and rerun power flow. In `run_fault_study.py`, `evaluate()` replays hourly schedules with unbalanced AC power flow at 15-minute intervals. Inspect `Fault_Results.csv` for the base-policy comparisons, `Severe_Stress_Results.csv` for longer outages, and `Severe_Tradeoff_Results.csv` for the energy bounds at maximum urgent-care coverage. For a fresh solve of one case in a separate extracted copy, use:

```sh
python severe_stress.py --event apr2024 --case long_reserve --policy coordinated --seconds 60
```

A fresh solve can stop at its time limit. Read the incumbent, bound, and AC-feasibility status before comparing it with the archived result.

## 5. Rebuild Fig. 4 and verify its voltage results

1. From the repository root, run `python fig4_final/replay/verify_voltage.py`. It checks the eight archived 15-minute voltage traces, the undervoltage hours in the 18–27 h backfeed zooms, and the April hour 22.75 comparison. The CSV `hour` column is the replay clock; subtract 12 to obtain hours from fault inception. The threshold is 0.95 p.u. and each sample represents 0.25 h.
2. To rebuild the exact curves and markers shown in the paper, use the saved TikZ source. Change into `fig4_final` and run `pdflatex fig4.tex`. The resulting page contains the full figure. Run `pdfcrop fig4.pdf fig4-cropped.pdf` if you need a tightly cropped figure for inclusion in another document. The figure source embeds its plotting coordinates, so this build requires no solver.
3. To recompute the newly added schedules, download and extract [the study bundle](reproduction.zip), install its requirements, and run the following commands from its root with the virtual environment active. Replace `PATH_TO_REPO` with the path to this repository. Repeat the three commands for `apr2024`.

   ```sh
   python PATH_TO_REPO/fig4_final/replay/run_matched_cases.py dec2023 long_no_surge autonomous
   python PATH_TO_REPO/fig4_final/replay/run_matched_cases.py dec2023 long_reserve autonomous
   python PATH_TO_REPO/fig4_final/replay/run_matched_cases.py dec2023 long_no_surge coordinated
   ```

   These commands write to `matched_case_runs/` under the extracted study directory. The surge joint care-first runs are in the original bundle and replay with `python severe_stress.py --replay --case long_reserve --policy coordinated`. A fresh solve can select a different schedule among ties. The archived CSVs are the inputs for checking the published values. The JSON and CSV files are simulated schedules and AC replay output.

## 6. Trace the other result figures

In [`build_severe_results.py`](build_severe_results.py), use `current_heatmap()` and `deadline_curve()` to trace equipment loading and missed urgent deadlines in Fig. 5. The inputs are the quarter-hour AC traces and treatment-start files under `extended_case_data/`. The generated sources are `figures/severe_current_heatmap.tikz` and `figures/severe_urgent_delays.tikz`.

The extracted bundle includes `README.txt` with model details, data provenance files, and the feeder's upstream notice in `feeder/License.txt`.
