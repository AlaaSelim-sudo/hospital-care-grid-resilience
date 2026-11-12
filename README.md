# Hospital Care Continuity Under Fault-Constrained Distribution Restoration

Numerical reproduction code for the IEEE PES General Meeting paper by Alaa Selim, Yitong Liu, Katie Lichter, and Junbo Zhao.

1. Download and extract [the numerical reproduction bundle](reproduction.zip).
2. In the extracted directory, use Python 3.12 and run:

   ```sh
   python -m pip install -r requirements.txt
   python reproduce.py --skip-pdf
   ```

The bundle contains the complete numerical workflow, archived solutions, solver bounds, aggregate public-source inputs and provenance, and TikZ result figures. The script checks the inputs and replays the archived schedules with OpenDSS. Hospital equipment, care cohorts, faults and switching times are modeled. See the included README.txt for model details and CITATION.cff for citation information.

The feeder files retain their upstream notice in feeder/License.txt.
