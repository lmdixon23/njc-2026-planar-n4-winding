#!/usr/bin/env python3
"""Certify two distinct preimages by the interval Krawczyk operator.

This is the artifact behind the "two certified preimages" proposition.  The
degree-winding route counts preimages globally and depends on the companion
uniform boundary-loop lemma; this certificate exhibits two of them directly and
depends on nothing beyond the matrices, the target, and outward-rounded
interval arithmetic.

Method
------
For a box X with midpoint m, let C = DF(m)^{-1} and

    K(X) = m - C F(m) + (I - C DF(X)) (X - m),

where DF(X) encloses the Jacobian over X.  If K(X) is contained in the
interior of X then F - y has exactly one zero in X (Krawczyk / Moore).
Disjointness of the two boxes then gives two distinct solutions.

Arithmetic is python-flint Arb ball arithmetic; every enclosure is outward, so
a reported PASS is rigorous at any working precision.  Lower precision merely
widens enclosures and can only turn a PASS into a FAIL.

Usage:
    python code/certify_preimage_boxes.py
    python code/certify_preimage_boxes.py --json-out results/preimage_boxes.json
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]

# Canonical network, exact rationals.
A_ROWS = [["1", "0", "-8", "2"], ["0", "1", "-1", "1/8"]]
B_ROWS = [["1", "0"], ["0", "1"], ["-1/5", "5/4"], ["-1/2", "9/5"]]
TARGET = ["-11/2", "1/10"]

# Approximate roots.  Their accuracy is irrelevant to correctness: they are box
# centres, and the Krawczyk test either certifies the box or does not.
CENTRES = [
    ("-0.56521578477251531355641846986664167934827388733396856",
     "3.11351292075998142167005832877332001362710351917071155"),
    ("25.17547640562957044984407684817497959995047290701810793",
     "7.4516846575961513667091955131549257381366376313620627"),
]
RADIUS = "1e-10"


def rat(text: str) -> arb:
    """Exact rational as a ball."""
    q = F(text)
    return arb(q.numerator) / arb(q.denominator)


def sigmoid(t: arb) -> arb:
    return 1 / (1 + (-t).exp())


def dsigmoid(t: arb) -> arb:
    s = sigmoid(t)
    return s * (1 - s)


def evaluate(A, B, x):
    """F(x) as a 2-vector of balls."""
    z = [B[i][0] * x[0] + B[i][1] * x[1] for i in range(4)]
    s = [sigmoid(zi) for zi in z]
    return [sum((A[r][i] * s[i] for i in range(4)), arb(0)) for r in range(2)]


def jacobian(A, B, x):
    """DF(x) as a 2x2 nested list of balls."""
    z = [B[i][0] * x[0] + B[i][1] * x[1] for i in range(4)]
    d = [dsigmoid(zi) for zi in z]
    return [[sum((A[r][i] * d[i] * B[i][c] for i in range(4)), arb(0))
             for c in range(2)] for r in range(2)]


def inverse2(M):
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    if not (det > 0) and not (det < 0):
        raise ValueError("centre Jacobian determinant not certified nonzero")
    return [[M[1][1] / det, -M[0][1] / det],
            [-M[1][0] / det, M[0][0] / det]], det


def krawczyk(A, B, y, centre, radius):
    """Return (certified, diagnostics) for the box of given radius about centre."""
    m = [rat(centre[0]), rat(centre[1])]
    r = rat(radius)
    # X as balls covering [m - r, m + r]
    X = [arb(m[0].mid(), r.mid()) if False else m[i] + arb(0, float(radius))
         for i in range(2)]

    Fm = evaluate(A, B, m)
    resid = [Fm[0] - y[0], Fm[1] - y[1]]
    Jm = jacobian(A, B, m)
    C, det_m = inverse2(Jm)
    JX = jacobian(A, B, X)

    # I - C JX
    E = [[(arb(1) if i == j else arb(0)) - sum((C[i][k] * JX[k][j] for k in range(2)), arb(0))
          for j in range(2)] for i in range(2)]
    # X - m  (a centred ball of radius r)
    d = [X[i] - m[i] for i in range(2)]

    K = []
    for i in range(2):
        acc = m[i]
        acc = acc - sum((C[i][k] * resid[k] for k in range(2)), arb(0))
        acc = acc + sum((E[i][k] * d[k] for k in range(2)), arb(0))
        K.append(acc)

    inside = []
    for i in range(2):
        lo_k, hi_k = K[i].lower(), K[i].upper()
        lo_x, hi_x = X[i].lower(), X[i].upper()
        inside.append(bool(lo_k > lo_x) and bool(hi_k < hi_x))

    return all(inside), {
        "centre": [str(m[0].str(30)), str(m[1].str(30))],
        "radius": radius,
        "residual_abs": [str(abs(resid[0]).upper()), str(abs(resid[1]).upper())],
        "centre_jacobian_det": str(det_m.str(20)),
        "krawczyk_image": [str(K[0].str(30, radius=True)), str(K[1].str(30, radius=True))],
        "strict_interior_by_coordinate": inside,
        "certified_unique_root": all(inside),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=256)
    parser.add_argument("--json-out")
    args = parser.parse_args()

    ctx.prec = args.precision

    A = [[rat(v) for v in row] for row in A_ROWS]
    B = [[rat(v) for v in row] for row in B_ROWS]
    y = [rat(v) for v in TARGET]

    certificates = []
    ok = True
    for centre in CENTRES:
        good, info = krawczyk(A, B, y, centre, RADIUS)
        certificates.append(info)
        ok = ok and good
        print(f"box at ({centre[0][:18]}, {centre[1][:18]}) "
              f"-> {'CERTIFIED' if good else 'NOT CERTIFIED'}")
        print(f"   residual upper bounds : {info['residual_abs']}")
        print(f"   centre det DF         : {info['centre_jacobian_det']}")
        print(f"   K(X) strictly inside X: {info['strict_interior_by_coordinate']}")

    # Disjointness: separation in the first coordinate.
    c0 = rat(CENTRES[0][0])
    c1 = rat(CENTRES[1][0])
    gap = c1 - c0 - 2 * rat(RADIUS)
    disjoint = bool(gap > 0)
    print(f"boxes disjoint: {disjoint}  (first-coordinate gap lower bound "
          f"{gap.lower()})")

    ok = ok and disjoint
    report = {
        "schema": "krawczyk-preimage-boxes-v1",
        "method": "interval Krawczyk operator, Arb ball arithmetic",
        "precision_bits": args.precision,
        "target_y": TARGET,
        "certificates": certificates,
        "boxes_disjoint": disjoint,
        "verdict": "PASS" if ok else "FAIL",
    }

    if args.json_out:
        out = ROOT / args.json_out
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2) + "\n",
                       encoding="utf-8", newline="\n")

    print(f"VERDICT: {'TWO DISTINCT PREIMAGES CERTIFIED' if ok else 'CERTIFICATION FAILED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
