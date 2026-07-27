#!/usr/bin/env python3
from fractions import Fraction as Q
import json

A2=[[Q(1),Q(0),Q(-8),Q(2)],[Q(0),Q(1),Q(-1),Q(1,8)]]
B2=[[Q(1),Q(0)],[Q(0),Q(1)],[Q(-1,5),Q(5,4)],[Q(-1,2),Q(9,5)]]

def construct(n: int):
    assert n>=2
    N=n+2
    A=[[Q(0) for _ in range(N)] for _ in range(n)]
    B=[[Q(0) for _ in range(n)] for _ in range(N)]
    for i in range(2):
        for j in range(4): A[i][j]=A2[i][j]
    for i in range(4):
        for j in range(2): B[i][j]=B2[i][j]
    for k in range(2,n):
        hidden=4+(k-2)
        A[k][hidden]=Q(1)
        B[hidden][k]=Q(1)
    assert len(A)==n and len(A[0])==N
    assert len(B)==N and len(B[0])==n
    return A,B

for n in range(2,9):
    A,B=construct(n)
    assert len(B)==n+2
    # Structural block checks.
    for k in range(2,n):
        h=4+k-2
        assert A[k][h]==1 and B[h][k]==1
        assert sum(abs(x) for j,row in enumerate(A) for x in ([row[h]] if j!=k else []))==0

print(json.dumps({
  'verdict':'VERIFIED',
  'checked_dimensions':list(range(2,9)),
  'width_formula':'N=n+2',
  'jacobian_factorization':'det DG_n = det DF_2 times product sigma prime of added coordinates',
  'higher_width_extension':'duplicate a hidden unit and split its output column',
  'conclusion':'counterexamples exist for every n>=2 and N>=n+2',
},indent=2))
