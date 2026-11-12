"""Replay archived schedules, validate the results, and rebuild the papers."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def run(*args):
    print('Running', *args, flush=True)
    subprocess.run(args, cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tectonic', default='tectonic', help='Path to the Tectonic executable')
    parser.add_argument('--skip-pdf', action='store_true', help='Validate results without compiling TeX')
    args = parser.parse_args()

    def python(script, *flags):
        run(sys.executable, script, *flags)

    python('analyze_evidence.py')
    python('run_fault_study.py')
    python('run_fault_study.py', '--sensitivity')
    python('fault_snapshots.py')
    python('annotate_results.py')
    python('audit_restoration_limits.py')
    python('build_fault_figures.py')
    python('validate_fault_study.py')
    python('verify_additional_results.py')
    python('compound_stress.py', '--replay')
    python('essential_reserve.py', '--replay')
    python('severe_stress.py', '--replay')
    python('severe_stress.py', '--reserve', '--replay')
    python('severe_stress.py', '--certificate', '--replay')
    if args.skip_pdf:
        python('build_compound_results.py')
        run(sys.executable, '-c', 'from build_severe_results import build; build()')
    else:
        python('build_fault_paper.py')
        compiler = shutil.which(args.tectonic) or str(Path(args.tectonic).resolve())
        run(compiler, '--keep-logs', 'Conference_Load_Figure.tex')
        run(compiler, '--keep-logs', 'main.tex')
        run(compiler, '--keep-logs', 'Fault_Study_Supplement.tex')
        shutil.copy2(ROOT / 'main.pdf', ROOT / 'PESGM_Fault_Care_Paper.pdf')


if __name__ == '__main__':
    main()
