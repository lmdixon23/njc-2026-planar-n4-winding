#!/usr/bin/env python3
from __future__ import annotations

import json
from functools import cmp_to_key
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parents[1] / 'results' / 'independent_sympy_reconstruction.json'
Q = sp.Rational

A = sp.Matrix([[1, 0, -8, 2], [0, 1, -1, Q(1, 8)]])
B = sp.Matrix([[1, 0], [0, 1], [Q(-1, 5), Q(5, 4)], [Q(-1, 2), Q(9, 5)]])
y = sp.Matrix([Q(-11, 2), Q(1, 10)])


def det2(u: sp.Matrix, v: sp.Matrix):
    return sp.expand(u[0] * v[1] - u[1] * v[0])


def half(v: sp.Matrix) -> int:
    return 0 if (v[1] > 0 or (v[1] == 0 and v[0] >= 0)) else 1


def angle_cmp(a, b):
    ha, hb = half(a[0]), half(b[0])
    if ha != hb:
        return -1 if ha < hb else 1
    c = det2(a[0], b[0])
    if c != 0:
        return -1 if c > 0 else 1
    return 0


def on_segment(q, a, b):
    return det2(b - a, q - a) == 0 and all(min(a[k], b[k]) <= q[k] <= max(a[k], b[k]) for k in range(2))


def winding(vertices, q):
    w = 0
    for i, a in enumerate(vertices):
        b = vertices[(i + 1) % len(vertices)]
        assert not on_segment(q, a, b)
        orient = det2(b - a, q - a)
        if a[1] <= q[1] < b[1] and orient > 0:
            w += 1
        elif b[1] <= q[1] < a[1] and orient < 0:
            w -= 1
    return w


def main():
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    out_minors = {}
    row_minors = {}
    mixed = {}
    for i, j in pairs:
        da = A[:, [i, j]].det()
        db = B[[i, j], :].det()
        out_minors[f'{i}{j}'] = da
        row_minors[f'{i}{j}'] = db
        mixed[f'{i}{j}'] = sp.factor(da * db)
    assert list(mixed.values()) == [1, Q(-5,4), Q(9,40), Q(8,5), -1, Q(53,200)]

    # Reconstruct the two row relations independently by solving linear systems.
    c1 = sp.Matrix.hstack(B[0, :].T, B[2, :].T).LUsolve(B[1, :].T)
    c2 = sp.Matrix.hstack(B[1, :].T, B[3, :].T).LUsolve(B[2, :].T)
    assert tuple(c1) == (Q(4,25), Q(4,5))
    assert tuple(c2) == (Q(53,100), Q(2,5))

    assert 15**1 * 5**4 * 4**5 > 8**1 * 4**4 * 5**5  # fifth-power form of group 1
    assert 245**47 > 94**47 * 2**53

    # Construct all positive and negative wall crossings from B, sort them exactly,
    # and recover the actual saturation polygon from a rational initial direction.
    walls = []
    for i in range(4):
        bx, by = B[i, 0], B[i, 1]
        plus = sp.Matrix([by, -bx])
        walls.append((plus, i, 1))
        walls.append((-plus, i, -1))
    walls = sorted(walls, key=cmp_to_key(angle_cmp))

    # Rotate cyclic order to start just after the rational direction d=(1,1/1000).
    d = sp.Matrix([1, Q(1,1000)])
    import math
    def angle_float(v):
        a = math.atan2(float(v[1]), float(v[0]))
        d0 = math.atan2(float(d[1]), float(d[0]))
        return (a - d0) % (2 * math.pi)
    walls = sorted(walls, key=lambda item: angle_float(item[0]))

    active = [i for i in range(4) if (B[i, :] * d)[0] > 0]
    current = sum((A[:, i] for i in active), sp.zeros(2, 1))
    vertices = []
    edge_sequence = []
    for _, i, sign in walls:
        vertices.append(current)
        edge = sign * A[:, i]
        edge_sequence.append((sign, i))
        current = current + edge
    assert current == vertices[0]
    w = winding(vertices, y)
    assert w == 2

    report = {
        'verdict': 'VERIFIED BY INDEPENDENT SYMPY RECONSTRUCTION',
        'mixed_coefficients': {k: str(v) for k, v in mixed.items()},
        'row_relations': [list(map(str, c1)), list(map(str, c2))],
        'edge_sequence': [[s, i] for s, i in edge_sequence],
        'vertices': [[str(p[0]), str(p[1])] for p in vertices],
        'target': [str(y[0]), str(y[1])],
        'winding': w,
        'scope': 'Independent internal reconstruction using SymPy matrices and a separately written wall-order and winding implementation. This is not external specialist review.',
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8', newline="\n")
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
