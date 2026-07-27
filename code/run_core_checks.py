#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CHECKS=[
    'code/verify_counterexample_exact.py',
    'code/two_negative_age_engine.py',
    'code/tropical_exposure_classifier.py',
    'code/check_alternating_classification_certificate.py',
    'code/verify_counterexample_family.py',
    'code/verify_dimension_lift.py',
    'code/verify_broad_analytic_family.py',
    'code/independent_reconstruction_sympy.py',
    'code/find_counterexample_preimages.py',
]

def main() -> None:
    for rel in CHECKS:
        path=ROOT/rel
        print(f'RUN {rel}',flush=True)
        subprocess.run([sys.executable,str(path)],cwd=ROOT,check=True)
    print('VERDICT: ALL CORE CHECKS PASSED')

if __name__=='__main__':
    main()
