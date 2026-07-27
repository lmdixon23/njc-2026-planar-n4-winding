#!/usr/bin/env python3
from fractions import Fraction
from itertools import combinations

G = [
    (Fraction(2), Fraction(3)),
    (Fraction(-4), Fraction(-4)),
    (Fraction(5), Fraction(-3)),
    (Fraction(-1), Fraction(2)),
]

EXPECTED = {
    (0, 1): Fraction(4),
    (0, 2): Fraction(-21),
    (0, 3): Fraction(7),
    (1, 2): Fraction(32),
    (1, 3): Fraction(-12),
    (2, 3): Fraction(7),
}

def det(u, v):
    return u[0] * v[1] - u[1] * v[0]

def main():
    minors = {(i, j): det(G[i], G[j]) for i, j in combinations(range(4), 2)}
    assert minors == EXPECTED

    d01 = minors[(0, 1)]
    d02 = minors[(0, 2)]
    d03 = minors[(0, 3)]
    d12 = minors[(1, 2)]
    d13 = minors[(1, 3)]
    d23 = minors[(2, 3)]

    # Gr(2,4) Plucker relation in the chosen sign convention.
    assert d01 * d23 - d02 * d13 + d03 * d12 == 0
    assert (-d02) * (-d13) == d01 * d23 + d03 * d12

    print('VERIFIED')
    print('output_minors =', minors)
    print('alternating_positive_magnitudes =', {
        'a': d01,
        'b': -d02,
        'c': d03,
        'd': d12,
        'e': -d13,
        'f': d23,
    })
    print('plucker =', (-d02) * (-d13), '=', d01 * d23 + d03 * d12)

if __name__ == '__main__':
    main()
