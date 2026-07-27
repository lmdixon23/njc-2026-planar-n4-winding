#!/usr/bin/env python3
from fractions import Fraction as Q
from itertools import combinations
import json

PAIRS=list(combinations(range(4),2))

def det(u,v): return u[0]*v[1]-u[1]*v[0]
def add(u,v): return (u[0]+v[0],u[1]+v[1])
def sub(u,v): return (u[0]-v[0],u[1]-v[1])
def scale(c,u): return (c*u[0],c*u[1])
def orient(a,b,c): return det(sub(b,a),sub(c,a))
def on_segment(x,a,b):
    return orient(a,b,x)==0 and min(a[0],b[0])<=x[0]<=max(a[0],b[0]) and min(a[1],b[1])<=x[1]<=max(a[1],b[1])

def proper_intersection_parameters(a,b,c,d):
    r=sub(b,a); ss=sub(d,c); den=det(r,ss)
    if den==0: return None
    ca=sub(c,a); t=det(ca,ss)/den; u=det(ca,r)/den
    return (t,u) if 0<t<1 and 0<u<1 else None

def exact_winding_spectrum(poly):
    split=[{Q(0),Q(1)} for _ in poly]
    for i,j in combinations(range(len(poly)),2):
        if j==(i+1)%len(poly) or i==(j+1)%len(poly):
            continue
        z=proper_intersection_parameters(poly[i],poly[(i+1)%len(poly)],poly[j],poly[(j+1)%len(poly)])
        if z:
            split[i].add(z[0]); split[j].add(z[1])
    M=max(abs(x) for p in poly for x in p)+10
    samples=[(M,M)]
    for i,a in enumerate(poly):
        b=poly[(i+1)%len(poly)]; e=sub(b,a); n=(-e[1],e[0])
        nn=max(abs(n[0]),abs(n[1]),Q(1))
        cuts=sorted(split[i])
        for lo,hi in zip(cuts,cuts[1:]):
            mid=add(a,scale((lo+hi)/2,e))
            for k in range(4,50):
                eps=Q(1,2**k)/nn
                q1=add(mid,scale(eps,n)); q2=sub(mid,scale(eps,n))
                if not any(on_segment(q1,poly[j],poly[(j+1)%len(poly)]) for j in range(len(poly))) and not any(on_segment(q2,poly[j],poly[(j+1)%len(poly)]) for j in range(len(poly))):
                    samples.extend([q1,q2]); break
            else:
                raise AssertionError('could not construct exact face samples')
    return sorted({winding(poly,q) for q in samples})

def winding(poly,q):
    w=0
    for i,a in enumerate(poly):
        b=poly[(i+1)%len(poly)]
        if on_segment(q,a,b): raise AssertionError('target lies on polygon')
        if a[1]<=q[1]:
            if b[1]>q[1] and orient(a,b,q)>0: w+=1
        else:
            if b[1]<=q[1] and orient(a,b,q)<0: w-=1
    return w

A=[(Q(1),Q(0)),(Q(0),Q(1)),(Q(-8),Q(-1)),(Q(2),Q(1,8))]
B=[(Q(1),Q(0)),(Q(0),Q(1)),(Q(-1,5),Q(5,4)),(Q(-1,2),Q(9,5))]
c=[Q(0)]*4

D={(i,j):det(A[i],A[j]) for i,j in PAIRS}
R={(i,j):det(B[i],B[j]) for i,j in PAIRS}
M={(i,j):D[(i,j)]*R[(i,j)] for i,j in PAIRS}

expected_D={(0,1):Q(1),(0,2):Q(-1),(0,3):Q(1,8),(1,2):Q(8),(1,3):Q(-2),(2,3):Q(1)}
expected_R={(0,1):Q(1),(0,2):Q(5,4),(0,3):Q(9,5),(1,2):Q(1,5),(1,3):Q(1,2),(2,3):Q(53,200)}
expected_M={(0,1):Q(1),(0,2):Q(-5,4),(0,3):Q(9,40),(1,2):Q(8,5),(1,3):Q(-1),(2,3):Q(53,200)}
assert D==expected_D
assert R==expected_R
assert M==expected_M
assert all(R[p]>0 for p in PAIRS)
assert D[(0,1)]*D[(2,3)]-D[(0,2)]*D[(1,3)]+D[(0,3)]*D[(1,2)]==0

# Exact logit relations used by the two AGE certificates.
# z1 = 4/25 z0 + 4/5 z2
assert add(scale(Q(4,25),B[0]),scale(Q(4,5),B[2]))==B[1]
# z2 = 53/100 z1 + 2/5 z3
assert add(scale(Q(53,100),B[1]),scale(Q(2,5),B[3]))==B[2]

# Split the shared positive coefficient m12 = 3/8 + 49/40.
assert Q(3,8)+Q(49,40)==M[(1,2)]

# Group 1 AM-GM coefficient:
# (5/4)^(4/5) (15/8)^(1/5) > 5/4 because 15/8 > 5/4.
assert Q(15,8)>Q(5,4)

# Group 2 AM-GM coefficient:
# (245/94)^(47/100) (1/2)^(53/100) > 1.
# Raise to the 100th power and check the exact integer inequality.
assert 245**47 > 94**47 * 2**53

# Actual row-wall traversal over one semicircle is +a1,+a2,+a3,-a0.
# Wall slopes are 0, p/q=4/25, r/s=5/18, infinity.
assert Q(0)<Q(4,25)<Q(5,18)
edges=[A[1],A[2],A[3],scale(Q(-1),A[0]),scale(Q(-1),A[1]),scale(Q(-1),A[2]),scale(Q(-1),A[3]),A[0]]
poly=[]
cur=A[0]
for e in edges:
    poly.append(cur)
    cur=add(cur,e)
assert cur==A[0]
expected_poly=[(Q(1),Q(0)),(Q(1),Q(1)),(Q(-7),Q(0)),(Q(-5),Q(1,8)),(Q(-6),Q(1,8)),(Q(-6),Q(-7,8)),(Q(2),Q(1,8)),(Q(0),Q(0))]
assert poly==expected_poly

y=(Q(-11,2),Q(1,10))
assert winding(poly,y)==2
spectrum=exact_winding_spectrum(poly)
assert spectrum==[0,1,2]

report={
  'verdict':'VERIFIED',
  'A_columns':[[str(x),str(y)] for x,y in A],
  'B_rows':[[str(x),str(y)] for x,y in B],
  'mixed_coefficients':{f'{i}{j}':str(M[(i,j)]) for i,j in PAIRS},
  'plucker':'verified',
  'age_group_1':'strictly dominates 5/4',
  'age_group_2_integer_check':f'{245**47} > {94**47 * 2**53}',
  'asymptotic_vertices':[[str(x),str(y)] for x,y in poly],
  'winding_target':[str(y[0]),str(y[1])],
  'winding':2,
  'winding_spectrum':spectrum,
  'conclusion':'For logistic and the stated activation class, det DF is everywhere positive while the degree-winding formula gives exactly two preimages of y.'
}
print(json.dumps(report,indent=2))
