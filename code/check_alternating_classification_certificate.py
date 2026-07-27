#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import Counter, deque
from fractions import Fraction
from pathlib import Path

import exact_alternating_classification as core

REPORT=Path(__file__).resolve().parents[1]/'results'/'exact_alternating_classification.json'

def parse_fraction(text: str) -> Fraction:
    if '/' in text:
        a,b=text.split('/',1)
        return Fraction(int(a),int(b))
    return Fraction(int(text),1)

def main() -> None:
    data=json.loads(REPORT.read_text(encoding='utf-8'))
    assert data['schema']=='exact-alternating-classification-v1'
    assert data['pattern_count']==48
    assert data['total_chambers']==2448
    assert data['pair_order']==['01','02','03','12','13','23']

    qualifying=[]
    total=0
    for rec in data['patterns']:
        D=[tuple(map(int,v)) for v in rec['directions']]
        pat=tuple(1 if core.det2(D[i],D[j])>0 else -1 for i,j in core.PAIRS)
        assert list(pat)==rec['pattern']
        W=core.walls_for(D)
        assert [list(w) for w in W]==rec['walls']
        normals=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]+W
        circuits=core.circuits(normals)
        assert len(circuits)==rec['circuit_count']

        chambers={}
        spectra=Counter()
        for item in rec['chambers']:
            signs=tuple(map(int,item['signs']))
            lam=tuple(parse_fraction(x) for x in item['lambda'])
            assert len(signs)==len(W)
            assert all(x>0 for x in lam)
            assert all(s*core.evalf(w,lam)>0 for s,w in zip(signs,W))
            assert not core.infeas(signs,circuits)
            chambers[signs]=lam
            spectra[core.spectrum(D,lam)]+=1

        assert len(chambers)==rec['chamber_count']
        total+=len(chambers)

        # Exact adjacency closure. Every one-wall flip is either another chamber
        # or is blocked by an exact positive circuit under Gordan's theorem.
        for signs in chambers:
            for i in range(len(W)):
                nbr=list(signs)
                nbr[i]*=-1
                nbr=tuple(nbr)
                assert nbr in chambers or core.infeas(nbr,circuits)

        # The feasible chamber graph must be connected.
        if chambers:
            start=next(iter(chambers))
            seen={start}
            queue=deque([start])
            while queue:
                signs=queue.popleft()
                for i in range(len(W)):
                    nbr=list(signs)
                    nbr[i]*=-1
                    nbr=tuple(nbr)
                    if nbr in chambers and nbr not in seen:
                        seen.add(nbr)
                        queue.append(nbr)
            assert len(seen)==len(chambers)

        encoded={'/'.join(map(str,k)):v for k,v in sorted(spectra.items())}
        assert encoded==rec['spectra']
        q=sum(n for spec,n in spectra.items() if max(spec)>=2 and min(spec)>=0)
        assert q==rec['qualifying_chambers']
        if q:
            qualifying.append(list(pat))

    assert total==2448
    assert qualifying==[[1,-1,1,1,-1,1]]
    assert data['qualifying_patterns']==qualifying
    print('VERDICT: PASS')
    print('patterns: 48')
    print('chambers: 2448')
    print('unique qualifying pattern: + - + + - +')

if __name__=='__main__':
    main()
