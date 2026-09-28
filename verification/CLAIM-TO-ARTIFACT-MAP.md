# Claim to artifact map

| Claim | Primary evidence | Separate check or failure condition |
|---|---|---|
| Exact minors, row circuits, and bridge inequalities | `code/verify_counterexample_exact.py` | Exact rational arithmetic; any sign or inequality failure exits nonzero |
| Everywhere-positive determinant | Two-bridge AGE argument, `code/two_negative_age_engine.py` | `code/tropical_exposure_classifier.py` separately classifies the exposure pattern |
| Asymptotic octagon and winding two | Exact integer orientation and ray-crossing predicates | `code/independent_reconstruction_sympy.py` rebuilds matrices, wall order, and winding in a separate internal implementation |
| Two distinct preimages / noninjectivity | `code/certify_preimage_boxes.py`, interval Krawczyk in Arb | Fixed nonsingular point preconditioner; strict inclusion in the exact rational boxes; boxes must be disjoint |
| Exact global multiplicity two | Winding two plus the companion degreeâ€“winding theorem | Requires the exact reviewed companion theorem and its hypotheses; not inferred from the two boxes alone |
| Explicit multiparameter counterexample family | `code/verify_broad_analytic_family.py` | Exact endpoint/monotonicity inequalities and winding-chamber conditions; `code/verify_counterexample_family.py` checks a one-parameter zero-bias slice |
| Dimension and width lift | `code/verify_dimension_lift.py` | Block determinant factorization checked symbolically; duplication preserves the represented function |
| Sharp threshold | Higher-width construction plus companion closed-width theorem at `N=n,n+1` | Companion theorem version must match the reviewed dependency |
| Symmetric log-concave activation class | Lemma 2.1 and the determinant decomposition in the manuscript | `code/two_negative_age_engine.py` checks the generic bridge algebra and coefficient inequalities; the infinite activation class is established analytically, not by finite sampling; exact multiplicity uses companion degreeâ€“winding |
| Canonical-realization chamber classification | `code/check_alternating_classification_certificate.py` | Complete unique canonical-pattern inventory, exact Gordan-circuit alternatives and adjacent-face sampling with exact clearance; theorem is restricted to open refined chambers of the enumerated canonical direction realizations |
| Repository integrity | `code/verify_sha256_manifest.py` | Rejects missing, extra, duplicated, malformed, or altered tracked paths |

## Scope

The conjecture refutation rests on the exact determinant argument together with
the two certified preimage boxes. It does **not** require the companion
degreeâ€“winding theorem. The companion paper is required for the stronger exact
global multiplicity statement and, through its `N=n` and `N=n+1` injectivity
theorem, for sharpness. The chamber classification is computer assisted and is
stated only on the open refined magnitude chambers of the canonical direction
realizations enumerated by its certificate. Wall magnitudes and unenumerated
direction realizations are outside that theorem.

## Method separation

Reports are recomputed rather than trusted: the manifest verifier reconstructs
every hash and compares the tracked inventory, and the preimage certificate is
evaluated in outward-rounded Arb ball arithmetic with a fixed point
preconditioner and exact rational box endpoints. The SymPy reconstruction is a
separate internal implementation, not external specialist review.
