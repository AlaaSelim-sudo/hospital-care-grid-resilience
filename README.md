# Hospital care and grid resilience under feeder faults

This study schedules emergency care and backup power together during a feeder outage and staged restoration. It follows two storm scenarios, December 2023 and April 2024, on a modified IEEE 13-node feeder with three modeled hospitals.

## Reproduce every figure in the submission

Clone this repository or download its source ZIP. Install a TeX compiler with IEEEtran, Latin Modern, TikZ, and PGFPlots. In a terminal at the repository root, extract [`reproduction.zip`](reproduction.zip) into a sibling folder named `pesgm-replay`. The two directories should be siblings and the repository directory should be named `hospital-care-grid-resilience`.

| Figure | Source used in the submission | How to reproduce it |
| --- | --- | --- |
| 1. Study workflow | `figures/care_grid_framework.tikz` in the study bundle | Render the archived hand-authored TikZ. |
| 2. NH/VT/ME county outage map | `figures/county_outage_maps.tikz` in the study bundle | Regenerate the exact TikZ from the archived county totals and boundaries, or render the saved source. |
| 3. Modified IEEE 13-node feeder | `figures/ieee13_single.tikz` in the study bundle | Render the archived hand-authored TikZ. |
| 4. Load, current, and voltage | [`fig4_final/figures/surge_replay.tikz`](fig4_final/figures/surge_replay.tikz) | Render the final TikZ and check the matched voltage traces as described below. |
| 5. Current heatmaps and care outcomes | `figures/severe_current_heatmap.tikz` and `figures/severe_urgent_delays.tikz` in the study bundle | Regenerate both sources from the archived simulations, or render the saved sources. |

### Render the figures

From the repository root, copy [`render_other_figures.tex`](figure_rebuild/render_other_figures.tex) into the extracted `pesgm-replay` folder. Run these commands in the specified directories. The tested compiler was Tectonic.

```sh
cd ../pesgm-replay
tectonic render_other_figures.tex
cd ../hospital-care-grid-resilience/fig4_final
tectonic fig4.tex
```

The first PDF contains Figs. 1, 2, 3, and 5 in that order. Its Fig. 5 wrapper sets `ConferenceResiliencePanels`, which selects the submission panels. The second PDF contains Fig. 4 on a full page. Use `pdfcrop fig4.pdf fig4-cropped.pdf` to trim its margins for inclusion elsewhere. The saved TikZ uses embedded plotting coordinates, so rendering does not need Python or an optimizer.

### Regenerate Fig. 2 and Fig. 5

From the repository root, use the following commands. For Fig. 5 use Python 3.12 with the extracted study's `requirements.txt` installed, as in Section 1.

```sh
python figure_rebuild/build_outage_map.py ../pesgm-replay
cd ../pesgm-replay
python ../hospital-care-grid-resilience/figure_rebuild/build_fig5.py
tectonic render_other_figures.tex
```

The map script combines the archived [county-event totals](figure_rebuild/County_Event_Power.csv) for the July 2023, December 2023, and April 2024 windows with `data/counties_NH_VT_ME.geojson` from the bundle. Totals are EAGLE-I recorded customer-outage-hours. Hatching indicates less than 75% coverage, and gray marks counties without retained records. The script starts from the archived county aggregates.

The Fig. 5 script reads the long-isolation 15-minute current and service replay files, archived treatment starts, and matching case summaries from the study bundle. It rebuilds the April equipment heatmaps, the December and April urgent deadline curves, and essential-service trajectories for 160 L and 320 L reserves. It uses archived schedules, so it does not rerun the optimizer.

The regenerated Fig. 2 and Fig. 5 sources were checked byte-for-byte against the submission. Figs. 1 and 3 come from the matching archived diagrams, and the rebuilt Fig. 4 page contents were checked against the submission figure. PDF page cropping and compiler metadata can differ.

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
2. To rebuild the exact curves and markers shown in the paper, use the saved TikZ source. Change into `fig4_final` and run `tectonic fig4.tex`. The resulting page contains the full figure. Run `pdfcrop fig4.pdf fig4-cropped.pdf` if you need a tightly cropped figure for inclusion in another document. The figure source embeds its plotting coordinates, so this build requires no solver.
3. To recompute the newly added schedules, download and extract [the study bundle](reproduction.zip), install its requirements, and run the following commands from its root with the virtual environment active. Replace `PATH_TO_REPO` with the path to this repository. Repeat the three commands for `apr2024`.

   ```sh
   python PATH_TO_REPO/fig4_final/replay/run_matched_cases.py dec2023 long_no_surge autonomous
   python PATH_TO_REPO/fig4_final/replay/run_matched_cases.py dec2023 long_reserve autonomous
   python PATH_TO_REPO/fig4_final/replay/run_matched_cases.py dec2023 long_no_surge coordinated
   ```

   These commands write to `matched_case_runs/` under the extracted study directory. The surge joint care-first runs are in the original bundle and replay with `python severe_stress.py --replay --case long_reserve --policy coordinated`. A fresh solve can select a different schedule among ties. The archived CSVs are the inputs for checking the published values. The JSON and CSV files are simulated schedules and AC replay output.

The extracted bundle includes `README.txt` with model details, data provenance files, and the feeder's upstream notice in `feeder/License.txt`.
