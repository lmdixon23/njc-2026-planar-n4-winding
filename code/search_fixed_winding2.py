#!/usr/bin/env python3
"""
Exploratory search for a positive normalized Jacobian margin with the exact
winding two output half edges from the source manuscript.

This script is not a proof. It uses floating point sampling and differential
evolution only to discover candidates or difficult negative witnesses.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import differential_evolution

G = np.array([
    [2.0, 3.0],
    [-4.0, -4.0],
    [5.0, -3.0],
    [-1.0, 2.0],
])

RESULT_PATH = Path(__file__).resolve().parents[1] / 'results' / 'fixed_winding2_initial.json'

def det(u: np.ndarray, v: np.ndarray) -> float:
    return float(u[0] * v[1] - u[1] * v[0])

def rho(z: np.ndarray) -> np.ndarray:
    a = np.abs(z)
    e = np.exp(-a)
    return e / (1.0 + e) ** 2

def make_points(seed: int = 0) -> np.ndarray:
    xs = np.linspace(-7.0, 7.0, 35)
    u, v = np.meshgrid(xs, xs, indexing='xy')
    grid = np.column_stack([u.ravel(), v.ravel()])
    rng = np.random.default_rng(seed)
    random_points = rng.uniform(-12.0, 12.0, size=(1200, 2))
    return np.vstack([grid, random_points])

POINTS = make_points()

def unpack(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    lp, lq, lr, ld, kappa, lam = theta
    p, q, r, delta = np.exp([lp, lq, lr, ld])
    s = (r * q + delta) / p

    B = np.array([
        [1.0, 0.0],
        [0.0, 1.0],
        [-p, q],
        [-r, s],
    ])
    c = np.array([0.0, 0.0, kappa, lam])
    values = {
        'p': float(p),
        'q': float(q),
        'r': float(r),
        's': float(s),
        'delta': float(delta),
        'kappa': float(kappa),
        'lambda': float(lam),
    }
    return B, c, values

def evaluate(theta: np.ndarray) -> dict[str, Any]:
    B, c, values = unpack(theta)
    z = POINTS @ B.T + c
    derivatives = rho(z)

    term_list = []
    pair_list = []
    mixed_list = []
    for i in range(4):
        for j in range(i + 1, 4):
            mixed = det(G[i], G[j]) * det(B[i], B[j])
            term = mixed * derivatives[:, i] * derivatives[:, j]
            term_list.append(term)
            pair_list.append([i, j])
            mixed_list.append(mixed)

    terms = np.asarray(term_list)
    jacobian = terms.sum(axis=0)
    denominator = np.abs(terms).sum(axis=0) + 1e-300
    ratio = jacobian / denominator
    index = int(np.argmin(ratio))

    return {
        'minimum_normalized_margin': float(ratio[index]),
        'minimum_jacobian_at_sample': float(jacobian[index]),
        'normalizing_sum_at_sample': float(denominator[index]),
        'witness_point': POINTS[index].tolist(),
        'row_parameters': values,
        'B': B.tolist(),
        'c': c.tolist(),
        'mixed_coefficients': [
            {'pair': pair, 'value': float(value)}
            for pair, value in zip(pair_list, mixed_list)
        ],
    }

def objective(theta: np.ndarray) -> float:
    return -evaluate(theta)['minimum_normalized_margin']

def main() -> None:
    bounds = [
        (-1.5, 1.5),
        (-1.5, 1.5),
        (-1.5, 1.5),
        (-2.5, 1.0),
        (-5.0, 5.0),
        (-5.0, 5.0),
    ]

    result = differential_evolution(
        objective,
        bounds=bounds,
        maxiter=25,
        popsize=8,
        tol=1e-5,
        polish=True,
        seed=2,
        workers=1,
    )

    report = evaluate(result.x)
    report.update({
        'status': 'candidate_not_found' if report['minimum_normalized_margin'] <= 0 else 'sampled_positive_candidate',
        'optimizer_success': bool(result.success),
        'optimizer_message': str(result.message),
        'optimizer_vector': result.x.tolist(),
        'warning': 'Floating point sampled evidence only. This is not a global certificate.',
    })

    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8', newline="\n")
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
