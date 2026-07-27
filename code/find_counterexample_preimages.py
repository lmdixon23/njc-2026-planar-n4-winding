#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
from scipy.optimize import root

A=np.array([[1.0,0.0,-8.0,2.0],[0.0,1.0,-1.0,0.125]])
B=np.array([[1.0,0.0],[0.0,1.0],[-0.2,1.25],[-0.5,1.8]])
y=np.array([-5.5,0.1])
OUT=Path(__file__).resolve().parents[1]/'results'/'counterexample_numerical_preimages.json'

def sigmoid(z: np.ndarray) -> np.ndarray:
    out=np.empty_like(z,dtype=float)
    pos=z>=0
    out[pos]=1.0/(1.0+np.exp(-z[pos]))
    ez=np.exp(z[~pos])
    out[~pos]=ez/(1.0+ez)
    return out

def F(x: np.ndarray) -> np.ndarray:
    return A@sigmoid(B@x)

def jacobian(x: np.ndarray) -> np.ndarray:
    s=sigmoid(B@x)
    return A@np.diag(s*(1-s))@B

def main() -> None:
    seeds=[]
    for u in np.linspace(-5,30,15):
        for v in np.linspace(-5,12,12):
            seeds.append(np.array([u,v]))
    roots=[]
    for seed in seeds:
        sol=root(lambda x:F(x)-y,seed,jac=lambda x:jacobian(x),method='hybr')
        if not sol.success or np.linalg.norm(F(sol.x)-y,np.inf)>1e-10:
            continue
        if not any(np.linalg.norm(sol.x-r,np.inf)<1e-7 for r in roots):
            roots.append(sol.x)
    roots.sort(key=lambda x:x[0])
    report={
      'status':'numerical corroboration only',
      'target':y.tolist(),
      'root_count_found':len(roots),
      'roots':[
        {
          'x':[float(z) for z in r],
          'F_of_x':[float(z) for z in F(r)],
          'infinity_residual':float(np.linalg.norm(F(r)-y,np.inf)),
          'det_DF':float(np.linalg.det(jacobian(r))),
        } for r in roots
      ],
      'proof_note':'The exact count two follows from the exact winding calculation and the degree-winding theorem, not from this numerical search.'
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8', newline="\n")
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
