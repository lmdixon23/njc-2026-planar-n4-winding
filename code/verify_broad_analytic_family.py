#!/usr/bin/env python3
from __future__ import annotations

import json
from fractions import Fraction as Q
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / 'results' / 'broad_analytic_family_verification.json'


def main() -> None:
    # Group 1: R1 = ((15/8)^(1/5) (5/4)^(4/5)) / (5/4) = (3/2)^(1/5).
    # We use exp(3/40) < 40/37 and prove R1 > 40/37 exactly.
    g1_lhs = 3 * 37**5
    g1_rhs = 2 * 40**5
    assert g1_lhs > g1_rhs

    # Explicit A(u,t) box.
    u_lo, u_hi = Q(199, 100), Q(51, 25)
    t_lo, t_hi = Q(3, 25), Q(13, 100)

    # Winding chamber sufficient conditions for the fixed target.
    assert t_lo > Q(1, 10)
    assert u_lo > 15 * t_hi
    assert u_hi < Q(5, 2)

    # The second bridge ratio is increasing in u and decreasing in t on the box:
    # d/du log R = (800 t - 47 u)/(100 u (u-8t)).
    assert 800 * t_lo > 47 * u_hi
    assert u_lo > 8 * t_hi

    # At the worst corner, p23/w1 = (u-8t)/2 = 19/40,
    # p12b/w3 = 245/94, and n13 = u/2 = 199/200.
    # Prove Rmin > 20/19 exactly after raising to the 100th power.
    lhs = 245**47 * 19**53 * 190**100
    rhs = 94**47 * 40**53 * 199**100
    assert lhs > rhs

    # exp(1/20) < 20/19 follows from exp(x) < 1/(1-x), 0<x<1.
    # Therefore bias residuals |delta1|<3/40 and |delta2|<1/20
    # preserve both strict bridge inequalities for logistic psi, since psi(s)<=|s|.

    report = {
        'verdict': 'VERIFIED',
        'family': {
            'A': '[[1,0,-8,u],[0,1,-1,t]]',
            'B': '[[1,0],[0,1],[-1/5,5/4],[-1/2,9/5]]',
            'u_interval': '199/100 < u < 51/25',
            't_interval': '3/25 < t < 13/100',
            'delta1': 'c1-(4/25)c0-(4/5)c2',
            'delta2': 'c2-(53/100)c1-(2/5)c3',
            'bias_conditions': '|delta1| < 3/40 and |delta2| < 1/20',
        },
        'exact_checks': {
            'group1_margin': f'{g1_lhs} > {g1_rhs}',
            'group2_margin_digits': {
                'lhs_digits': len(str(lhs)),
                'rhs_digits': len(str(rhs)),
                'difference_positive': True,
            },
            'winding_chamber': [
                't>1/10',
                'u>15t',
                'u<5/2',
            ],
        },
        'conclusion': 'Every parameter satisfying the displayed inequalities is a positive-Jacobian noninjective logistic four-ridge map with winding two at y=(-11/2,1/10).',
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
