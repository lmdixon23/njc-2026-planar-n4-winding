#!/usr/bin/env python3
"""Certify two distinct preimages by the interval Krawczyk operator.

This is the artifact behind the "two certified preimages" proposition.  The
degree-winding route counts preimages globally and depends on the companion
uniform boundary-loop lemma; this certificate exhibits two of them directly and
depends on nothing beyond the matrices, the target, and outward-rounded
interval arithmetic.

Method
------
For a box X with midpoint m, let C be a fixed point matrix approximating DF(m)^(-1) and

    G(x) = F(x) - y,
    K_G(X) = m - C G(m) + (I - C DG(X)) (X - m),

where DG(X)=DF(X) encloses the Jacobian over X.  If K_G(X) is contained
in the interior of X then G has exactly one zero in X (Krawczyk / Moore).
Disjointness of the two boxes then gives two distinct solutions.

Arithmetic is python-flint Arb ball arithmetic; every enclosure is outward, so
a reported PASS is rigorous at the working precision.  Insufficient precision
may fail to certify an otherwise valid box.

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
    ("-0.5652157847725153", "3.1135129207599814"),
    ("25.1754764056295704", "7.4516846575961514"),
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
    centre_q = [F(centre[0]), F(centre[1])]
    radius_q = F(radius)
    m = [rat(centre[0]), rat(centre[1])]

    # X encloses the exact rational box [centre-radius, centre+radius].
    # Arb rounds the radius outward; this can only widen the derivative enclosure.
    X = [m[i] + arb(0, radius) for i in range(2)]

    Fm = evaluate(A, B, m)
    resid = [Fm[0] - y[0], Fm[1] - y[1]]
    Jm = jacobian(A, B, m)
    C_ball, det_m = inverse2(Jm)
    # The Krawczyk theorem uses a fixed point preconditioner.  Take the dyadic
    # midpoint of the outward-rounded inverse enclosure, then certify it is
    # nonsingular.  This avoids treating an interval matrix as the theorem's C.
    C = [[arb(C_ball[i][j].mid()) for j in range(2)] for i in range(2)]
    det_C = C[0][0] * C[1][1] - C[0][1] * C[1][0]
    if not (det_C > 0) and not (det_C < 0):
        raise ValueError("fixed Krawczyk preconditioner not certified nonsingular")
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
    exact_box = []
    for i in range(2):
        lo_q = centre_q[i] - radius_q
        hi_q = centre_q[i] + radius_q
        lo_x = rat(f"{lo_q.numerator}/{lo_q.denominator}")
        hi_x = rat(f"{hi_q.numerator}/{hi_q.denominator}")
        exact_box.append([str(lo_q), str(hi_q)])
        lo_k, hi_k = K[i].lower(), K[i].upper()
        # Compare with the exact rational endpoints, allowing only certified
        # strict separation from their outward-rounded Arb enclosures.
        inside.append(bool(lo_k > lo_x.upper()) and bool(hi_k < hi_x.lower()))

    return all(inside), {
        "centre": [centre[0], centre[1]],
        "centre_ball": [str(m[0].str(30)), str(m[1].str(30))],
        "radius": radius,
        "exact_box_endpoints": exact_box,
        "preconditioner_det": str(det_C.str(20)),
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
    parser.add_argument("--check-json")
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
        "schema": "krawczyk-preimage-boxes-v2",
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

    if args.check_json:
        expected_path = ROOT / args.check_json
        expected = json.loads(expected_path.read_text(encoding="utf-8"))
        if expected != report:
            print(f"canonical report mismatch: {args.check_json}")
            return 1
        print(f"canonical report matches: {args.check_json}")

    print(f"VERDICT: {'TWO DISTINCT PREIMAGES CERTIFIED' if ok else 'CERTIFICATION FAILED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
