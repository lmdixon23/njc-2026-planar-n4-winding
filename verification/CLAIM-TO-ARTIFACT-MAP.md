# Claim to artifact map

| Claim | Primary evidence | Independent check or failure condition |
|---|---|---|
| Exact minors, row circuits, and bridge inequalities | `code/verify_counterexample_exact.py` | Exact rational arithmetic; any sign or inequality failure exits nonzero |
| Everywhere-positive determinant | Two-bridge AGE argument, `code/two_negative_age_engine.py` | `code/tropical_exposure_classifier.py` classifies the exposure pattern independently |
| Asymptotic octagon and winding two | Exact integer orientation and ray-crossing predicates | `code/independent_reconstruction_sympy.py` rebuilds matrices, wall order, and winding in a separate implementation |
| Two distinct preimages | `code/certify_preimage_boxes.py`, interval Krawczyk in Arb | Krawczyk image must lie strictly inside each box coordinatewise, and the boxes must be disjoint; independent of the degree-winding route |
| Multiparameter counterexample family | `code/verify_counterexample_family.py` | Exact rational comparisons over the parameter rectangle |
| Dimension and width lift | `code/verify_dimension_lift.py` | Block determinant factorization checked symbolically |
| Broad activation class | `code/verify_broad_analytic_family.py` | Bridge inequality re-derived for each admissible potential |
| Alternating chamber classification | `code/check_alternating_classification_certificate.py` | Exact Gordan-circuit alternatives; the theorem is stated over the enumerated canonical realizations only |
| Repository integrity | `code/verify_sha256_manifest.py` | Rejects missing, extra, duplicated, malformed, or altered tracked paths |

## Scope

The refutation rests on the exact determinant argument together with the
certified preimage boxes. The degree-winding computation supplies the exact
global multiplicity and explains the count, but is not required for
noninjectivity. The chamber classification is computer assisted and is stated
over the canonical direction realizations that its certificate enumerates.

## Independence

Reports are recomputed rather than trusted: the manifest verifier reconstructs
every hash and compares the tracked inventory, and the preimage certificate is
evaluated in outward-rounded ball arithmetic so that a reported PASS is
rigorous at any working precision.
