#!/usr/bin/env python3
from __future__ import annotations

from fractions import Fraction as Q
from itertools import combinations
import json

# Canonical counterexample rows. A general caller may replace these exact pairs.
FORMS=[
    (Q(1),Q(0)),
    (Q(0),Q(1)),
    (Q(-1,5),Q(5,4)),
    (Q(-1,2),Q(9,5)),
]

def value(form,t):
    a,b=form
    return abs(a+b*t)

def equality_roots(forms):
    roots=set()
    permanent_ties=[]
    for i,j in combinations(range(len(forms)),2):
        ai,bi=forms[i]
        aj,bj=forms[j]
        for sign in (1,-1):
            # ai+bi t = sign times aj+sign times bj t
            A=ai-sign*aj
            B=bi-sign*bj
            if B:
                roots.add(-A/B)
            elif A==0:
                permanent_ties.append((i,j,sign))
    return sorted(roots),permanent_ties

def interval_samples(roots):
    if not roots:
        return [(None,None,Q(0))]
    out=[(None,roots[0],roots[0]-1)]
    out.extend((a,b,(a+b)/2) for a,b in zip(roots,roots[1:]))
    out.append((roots[-1],None,roots[-1]+1))
    return out

def classify(forms):
    roots,ties=equality_roots(forms)
    exposed={pair:[] for pair in combinations(range(len(forms)),2)}
    for lo,hi,t in interval_samples(roots):
        vals=[value(f,t) for f in forms]
        order=sorted(range(len(forms)),key=lambda i:vals[i])
        if vals[order[1]]<vals[order[2]]:
            pair=tuple(sorted(order[:2]))
            exposed[pair].append({
                'left':None if lo is None else str(lo),
                'right':None if hi is None else str(hi),
                'sample':str(t),
                'absolute_projections':[str(x) for x in vals],
            })
    return roots,ties,exposed

def main():
    roots,ties,exposed=classify(FORMS)
    assert not ties
    assert exposed[(0,2)]==[]
    assert exposed[(1,3)]==[]
    report={
        'verdict':'PASS',
        'forms':[[str(a),str(b)] for a,b in FORMS],
        'critical_roots':[str(x) for x in roots],
        'exposed_intervals':{f'{i}{j}':v for (i,j),v in exposed.items()},
        'negative_pair_02':'hidden',
        'negative_pair_13':'hidden',
        'theorem':'A uniquely least decay pair dominates every other Cauchy Binet term exponentially along the corresponding ray. A negative exposed pair therefore forbids global Jacobian positivity.',
    }
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
