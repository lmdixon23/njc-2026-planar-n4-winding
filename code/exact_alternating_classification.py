from fractions import Fraction
from itertools import permutations, product, combinations
from collections import Counter
from functools import reduce
from math import gcd
import math, json, time, hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

BASE=[(1,0),(1,1),(0,1),(-1,1)]
PAIRS=list(combinations(range(4),2))

def det2(a,b): return a[0]*b[1]-a[1]*b[0]
patterns={}
for perm in permutations(range(4)):
    for flips in product([1,-1],repeat=4):
        D=[(flips[i]*BASE[perm[i]][0],flips[i]*BASE[perm[i]][1]) for i in range(4)]
        pat=tuple(1 if det2(D[i],D[j])>0 else -1 for i,j in PAIRS)
        patterns.setdefault(pat,D)

def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def scale(c,a):return tuple(c*x for x in a)
Z=(Fraction(0),)*4

def norm(coeff):
    den=1
    for c in coeff:den=math.lcm(den,c.denominator)
    ints=[int(c*den) for c in coeff]
    g=reduce(gcd,[abs(x) for x in ints])
    if g:ints=[x//g for x in ints]
    for x in ints:
        if x:
            if x<0:ints=[-y for y in ints]
            break
    return tuple(ints)

def walls_for(D):
    edges=[];h=[]
    for sg in [1,-1]:
        for i,d in enumerate(D):
            x=[Fraction(0)]*4;y=[Fraction(0)]*4
            x[i]=sg*d[0];y[i]=sg*d[1]
            edges.append((tuple(x),tuple(y)));h.append((sg*d[0],sg*d[1]))
    verts=[];cur=(Z,Z)
    for e in edges:
        verts.append(cur);cur=(add(cur[0],e[0]),add(cur[1],e[1]))
    L=set();Ms={};Ds={}
    for r in range(8):
        for s in range(8):
            if s in {r,(r+1)%8}:continue
            dx=sub(verts[s][0],verts[r][0]);dy=sub(verts[s][1],verts[r][1])
            c=sub(scale(Fraction(h[r][0]),dy),scale(Fraction(h[r][1]),dx))
            if any(c):L.add(norm(c))
    for r in range(8):
        for s in range(8):
            if r==s:continue
            drs=det2(h[r],h[s])
            if drs==0:continue
            dx=sub(verts[s][0],verts[r][0]);dy=sub(verts[s][1],verts[r][1])
            M=sub(scale(Fraction(h[s][1]),dx),scale(Fraction(h[s][0]),dy))
            Ms[(r,s)]=M;Ds[(r,s)]=Fraction(drs)
    K=set()
    for r in range(8):
        comps=[s for s in range(8) if (r,s) in Ms]
        for s,t in combinations(comps,2):
            c=sub(scale(Ds[(r,t)],Ms[(r,s)]),scale(Ds[(r,s)],Ms[(r,t)]))
            if any(c):K.add(norm(c))
    return [w for w in sorted(L|K) if min(w)<0<max(w)]

def det_bareiss(M):
    n=len(M)
    if n==0:return 1
    A=[list(map(int,row)) for row in M];sg=1;prev=1
    for k in range(n-1):
        if A[k][k]==0:
            sw=next((i for i in range(k+1,n) if A[i][k]!=0),None)
            if sw is None:return 0
            A[k],A[sw]=A[sw],A[k];sg=-sg
        piv=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):A[i][j]=(A[i][j]*piv-A[i][k]*A[k][j])//prev
        prev=piv
        for i in range(k+1,n):A[i][k]=0
    return sg*A[-1][-1]

def circuit_coeffs(rows):
    k=len(rows)
    for cols in combinations(range(4),k-1):
        S=[[r[j] for j in cols] for r in rows]
        c=[]
        for i in range(k):c.append(((-1)**i)*det_bareiss([S[j] for j in range(k) if j!=i]))
        if not all(c):continue
        if all(sum(c[i]*rows[i][j] for i in range(k))==0 for j in range(4)):
            g=reduce(gcd,[abs(x) for x in c]);c=[x//g for x in c]
            if c[0]<0:c=[-x for x in c]
            return tuple(c)
    return None

def circuits(normals):
    out=[]
    for k in range(2,min(5,len(normals))+1):
        for idx in combinations(range(len(normals)),k):
            c=circuit_coeffs([normals[i] for i in idx])
            if c:out.append((idx,c))
    return out

def evalf(w,v):return sum(Fraction(w[i])*v[i] for i in range(4))
def init_sign(W):
    for q in [(1,2,4,8),(1,3,9,27),(2,3,5,7),(1,5,2,11)]:
        v=tuple(map(Fraction,q))
        if all(evalf(w,v)!=0 for w in W):return tuple(1 if evalf(w,v)>0 else -1 for w in W)
    raise RuntimeError

def infeas(signs,C):
    S=(1,1,1,1)+tuple(signs)
    for idx,c in C:
        z=[c[k]*S[idx[k]] for k in range(len(idx))]
        if all(x>0 for x in z) or all(x<0 for x in z):return True
    return False

def sample(signs,W):
    rows=np.array([(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]+[tuple(s*x for x in w) for s,w in zip(signs,W)],float)
    A=np.column_stack([-rows,np.ones(len(rows))])
    R=linprog([0,0,0,0,-1],A_ub=A,b_ub=np.zeros(len(rows)),A_eq=[[1,1,1,1,0]],b_eq=[1],bounds=[(None,None)]*5,method='highs')
    if not R.success or R.x[4]<1e-10:raise RuntimeError('Gordan/LP disagreement')
    for dm in [10**3,10**5,10**7,10**9]:
        v=tuple(Fraction(float(x)).limit_denominator(dm) for x in R.x[:4]);ss=sum(v);v=tuple(x/ss for x in v)
        if all(x>0 for x in v) and all(s*evalf(w,v)>0 for s,w in zip(signs,W)):return v
    raise RuntimeError('rationalize')

def enumerate_chambers(W):
    N=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]+W;C=circuits(N)
    s0=init_sign(W);assert not infeas(s0,C);v0=sample(s0,W)
    Q=[s0];seen={s0:v0};blocked=set()
    while Q:
        s=Q.pop()
        for i in range(len(W)):
            z=list(s);z[i]*=-1;z=tuple(z)
            if z in seen or z in blocked:continue
            if infeas(z,C):blocked.add(z)
            else:seen[z]=sample(z,W);Q.append(z)
    # exact closure check
    assert len(seen)>0
    for s in seen:
        for i in range(len(W)):
            z=list(s);z[i]*=-1;z=tuple(z)
            assert z in seen or infeas(z,C)
    return seen,C

# exact planar winding

def p_add(a,b):return (a[0]+b[0],a[1]+b[1])
def p_sub(a,b):return (a[0]-b[0],a[1]-b[1])
def p_scale(c,a):return (c*a[0],c*a[1])
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def orient(a,b,c):return cross(p_sub(b,a),p_sub(c,a))
def on(q,a,b):return orient(a,b,q)==0 and min(a[0],b[0])<=q[0]<=max(a[0],b[0]) and min(a[1],b[1])<=q[1]<=max(a[1],b[1])
def ipar(a,b,c,d):
    r=p_sub(b,a);s=p_sub(d,c);den=cross(r,s)
    if den==0:return None
    ca=p_sub(c,a);t=cross(ca,s)/den;u=cross(ca,r)/den
    return (t,u) if 0<t<1 and 0<u<1 else None

def wind(P,q):
    w=0
    for i in range(8):
        a=P[i];b=P[(i+1)%8]
        if on(q,a,b):raise RuntimeError('sample on edge')
        if a[1]<=q[1]:
            if b[1]>q[1] and orient(a,b,q)>0:w+=1
        else:
            if b[1]<=q[1] and orient(a,b,q)<0:w-=1
    return w

def spectrum(D,lam):
    E=[(lam[i]*D[i][0],lam[i]*D[i][1]) for i in range(4)];E+= [(-x,-y) for x,y in E]
    P=[];cur=(Fraction(0),Fraction(0))
    for e in E:P.append(cur);cur=p_add(cur,e)
    split=[{Fraction(0),Fraction(1)} for _ in range(8)]
    for i,j in combinations(range(8),2):
        if j==(i+1)%8 or i==(j+1)%8:continue
        z=ipar(P[i],P[(i+1)%8],P[j],P[(j+1)%8])
        if z:split[i].add(z[0]);split[j].add(z[1])
    M=max(abs(x) for p in P for x in p)+10;samples=[(M,M)]
    for i in range(8):
        a=P[i];b=P[(i+1)%8];e=p_sub(b,a);n=(-e[1],e[0]);nn=max(abs(n[0]),abs(n[1]),Fraction(1))
        T=sorted(split[i])
        for x,y in zip(T,T[1:]):
            m=p_add(a,p_scale((x+y)/2,e))
            for k in range(4,40):
                ep=Fraction(1,2**k)/nn;q1=p_add(m,p_scale(ep,n));q2=p_sub(m,p_scale(ep,n))
                if not any(on(q1,P[j],P[(j+1)%8]) for j in range(8)) and not any(on(q2,P[j],P[(j+1)%8]) for j in range(8)):
                    samples.extend([q1,q2]);break
            else:raise RuntimeError('perturb')
    return tuple(sorted({wind(P,q) for q in samples}))

def fs(x):return f'{x.numerator}/{x.denominator}'

def main():
    OUT=[];qual=[];total=0;t0=time.time()
    for no,(pat,D) in enumerate(sorted(patterns.items()),1):
        W=walls_for(D);CH,C=enumerate_chambers(W);counts=Counter()
        for v in CH.values():counts[spectrum(D,v)]+=1
        Q=sum(n for s,n in counts.items() if max(s)>=2 and min(s)>=0)
        rec={'pattern':list(pat),'directions':[list(x) for x in D],'walls':[list(x) for x in W],
             'circuit_count':len(C),'chamber_count':len(CH),'spectra':{'/'.join(map(str,k)):v for k,v in sorted(counts.items())},
             'qualifying_chambers':Q,
             'chambers':[{'signs':list(s),'lambda':[fs(x) for x in v]} for s,v in sorted(CH.items())]}
        OUT.append(rec);total+=len(CH)
        if Q:qual.append(rec['pattern'])
        print(no,pat,len(W),len(CH),dict(counts),'Q',Q,flush=True)
    report={'schema':'exact-alternating-classification-v1','pattern_count':len(OUT),'total_chambers':total,
            'qualifying_patterns':qual,'result':'unique qualifying chirotope' if len(qual)==1 else 'unexpected',
            'pair_order':['01','02','03','12','13','23'],
            'method':'48 exact rank-two chirotopes; refined L and K magnitude walls; exact Gordan circuit completeness; rational chamber witnesses; exact rational winding spectra',
            'patterns':OUT,'elapsed_seconds':time.time()-t0}
    raw=json.dumps(report,indent=2).encode();out_path=Path(__file__).resolve().parents[1]/'results'/'exact_alternating_classification.json'
    out_path.parent.mkdir(parents=True,exist_ok=True)
    out_path.write_bytes(raw)
    print('DONE',len(OUT),total,qual,time.time()-t0,'sha',hashlib.sha256(raw).hexdigest())

if __name__ == '__main__':
    main()
