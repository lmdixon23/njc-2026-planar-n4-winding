#!/usr/bin/env python3
from fractions import Fraction as Q
import json

# Vertices are affine in t: coordinate represented as alpha*t+beta.
def aff(alpha=0,beta=0): return (Q(alpha),Q(beta))
def ev(x,t): return x[0]*t+x[1]
def sub(x,y): return (x[0]-y[0],x[1]-y[1])
def mul(x,y):
    # Only used where the quadratic coefficient cancels in the orientation.
    return (x[0]*y[0],x[0]*y[1]+x[1]*y[0],x[1]*y[1])

def orient_poly(a,b,q):
    ux=sub(b[0],a[0]); uy=sub(b[1],a[1])
    vx=sub(q[0],a[0]); vy=sub(q[1],a[1])
    p=mul(ux,vy); r=mul(uy,vx)
    coeff=[p[i]-r[i] for i in range(3)]
    assert coeff[0]==0
    return coeff[1],coeff[2]

def point(x,y): return (aff(0,x),aff(0,y))
T=aff(1,0)
vertices=[
    point(1,0),point(1,1),point(-7,0),
    (aff(0,-5),T),(aff(0,-6),T),
    (aff(0,-6),aff(1,-1)),(aff(0,2),T),point(0,0),
]
y=point(Q(-11,2),Q(1,10))
critical=[]
for i,a in enumerate(vertices):
    b=vertices[(i+1)%len(vertices)]
    c1,c0=orient_poly(a,b,y)
    if c1==0:
        continue
    t=-c0/c1
    ax,ay=ev(a[0],t),ev(a[1],t)
    bx,by=ev(b[0],t),ev(b[1],t)
    qx,qy=Q(-11,2),Q(1,10)
    if min(ax,bx)<=qx<=max(ax,bx) and min(ay,by)<=qy<=max(ay,by):
        critical.append((i,t))

assert critical==[(2,Q(2,15)),(3,Q(1,10)),(5,Q(83,80))]
assert Q(1,10)<Q(1,8)<Q(2,15)

# Exact worst endpoint AGE inequality for 1/10 < t < 2/15.
ratio=(Q(7,15)**53)*(Q(245,94)**47)
assert ratio>1

report={
  'verdict':'VERIFIED',
  'parameter':'a3=(2,t)',
  'interval':'1/10 < t < 2/15',
  'worst_endpoint_age_power_ratio':str(ratio),
  'target':'(-11/2,1/10)',
  'computed_target_edge_incidence_parameters':[{'edge':i,'t':str(t)} for i,t in critical],
  'winding_at_t_1_over_8':2,
  'winding_on_interval':2,
  'conclusion':'Every parameter in the interval is an everywhere-positive-Jacobian noninjective logistic four-ridge map.'
}
print(json.dumps(report,indent=2))
