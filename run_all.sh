#!/usr/bin/env bash
# Reproduce the maintained deterministic claims in mathematical order.
set -euo pipefail
export PYTHONHASHSEED=0

# Select an interpreter that carries the pinned dependencies, not merely one
# that starts. An explicit choice may be supplied as PY=/path/to/python.
# On systems where `bash` resolves to the Windows Subsystem for Linux, the
# default interpreter is a separate installation and may lack them.
PY="${PY:-python3}"
"$PY" -c "import flint" 2>/dev/null || PY=python
"$PY" -c "import flint" 2>/dev/null || {
  echo "no interpreter with the pinned dependencies found (python-flint missing)."
  echo "install requirements-lock.txt, or run: PY=/path/to/python bash run_all.sh"
  exit 1
}
"$PY" -c "import sys, flint; print('Interpreter:', sys.executable, '| flint', flint.__version__)"

printf '%s\n' '== Canonical-file integrity =='
"$PY" code/verify_sha256_manifest.py

printf '%s\n' '== Exact counterexample: minors, circuits, bridges =='
"$PY" code/verify_counterexample_exact.py

printf '%s\n' '== Two-negative AGE engine =='
"$PY" code/two_negative_age_engine.py

printf '%s\n' '== Tropical exposure classification =='
"$PY" code/tropical_exposure_classifier.py

printf '%s\n' '== Alternating chamber certificate =='
"$PY" code/test_classification_contracts.py
"$PY" code/check_alternating_classification_certificate.py

printf '%s\n' '== One-parameter counterexample slice =='
"$PY" code/verify_counterexample_family.py

printf '%s\n' '== Dimension and width lift =='
"$PY" code/verify_dimension_lift.py

printf '%s\n' '== Explicit multiparameter counterexample family =='
"$PY" code/verify_broad_analytic_family.py

printf '%s\n' '== Certified preimage boxes (interval Krawczyk) =='
"$PY" code/certify_preimage_boxes.py --check-json results/preimage_boxes.json

printf '%s\n' '== Separate internal SymPy reconstruction =='
"$PY" code/independent_reconstruction_sympy.py

printf '%s\n' '== Final canonical-file integrity =='
"$PY" code/verify_sha256_manifest.py

printf '%s\n' 'VERDICT: ALL MAINTAINED CHECKS PASSED'
