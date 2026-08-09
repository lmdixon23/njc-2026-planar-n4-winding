# A positive-Jacobian noninjective sigmoid ridge map and the sharp hidden-width threshold

![verify](https://github.com/lmdixon23/njc-2026-planar-n4-winding/actions/workflows/verify.yml/badge.svg)

Manuscript source, exact data, and reproducibility materials for an exact
planar four-ridge counterexample to the Neural Jacobian Conjecture and the
resulting sharp hidden-width threshold.

## Main result

At n = 2 and N = 4, the manuscript gives exact rational matrices defining a sigmoid ridge map with positive Jacobian determinant everywhere that is nevertheless noninjective. An explicit rational target has winding number two, and two distinct preimages are independently certified by interval Krawczyk methods. Together with the positive results at widths N = n and N = n + 1, this gives the sharp hidden-width threshold.

## Repository contents

- `paper/` LaTeX manuscript.
- `data/` the exact rational counterexample, byte frozen.
- `code/` maintained certificate and reconstruction scripts.
- `results/` machine-readable reports referenced by the manuscript.
- `verification/` claim-to-artifact map.

Working notes, exploratory searches, and superseded drafts live in the ignored
`_local/` directory and are not part of the public repository.

## Maintained verification

```
python -m pip install -r requirements-lock.txt
bash run_all.sh
```

The expected final line is

```
VERDICT: ALL MAINTAINED CHECKS PASSED
```

The gate must run under the interpreter carrying `requirements-lock.txt`. The
script selects one that can import `python-flint` and prints it before the
first check; an explicit choice may be supplied as
`PY=/path/to/python bash run_all.sh`. On systems where `bash` resolves to the
Windows Subsystem for Linux, the default interpreter is a separate
installation and may lack those dependencies.

## What the gate establishes

The determinant of the counterexample is positive everywhere by a two-bridge
arithmetic-geometric-exponential argument, checked in exact rational
arithmetic. Its asymptotic octagon has winding number two at an explicit
rational target, computed by exact integer predicates. Two distinct preimages
of that target are certified independently by the interval Krawczyk operator
in Arb ball arithmetic, so the refutation does not depend on the degree-winding
route. The multiparameter family, the dimension and width lift, the broad
activation class, and the alternating chamber classification are checked
separately.

## Canonical data

`data/counterexample.json` carries the exact rational matrices and target. Its
SHA-256 is recorded in `SHA256SUMS.txt` and bound into every report.

## Related research

The companion [asymptotic-polygon global-inversion repository](https://github.com/lmdixon23/njc-2026-cluster-inversion) develops the planar asymptotic-polygon and winding framework used in the surrounding Neural Jacobian programme.

## Author

Logan M. Dixon · [research site](https://lmdixon23.github.io/) · [ORCID](https://orcid.org/0009-0001-0592-462X)

## Build

```
cd paper
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```
