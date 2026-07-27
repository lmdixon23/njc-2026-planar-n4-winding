#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

B=np.array([[1.0,0.0],[0.0,1.0],[-0.2,1.25],[-0.5,1.8]])
COEFF=np.array([1.0,-1.25,0.225,1.6,-1.0,0.265])
PAIRS=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
OUT=Path(__file__).resolve().parents[1]/'results'/'numerical_stress_test.json'

def log_rho(z):
    a=np.abs(z)
    return -a-2*np.log1p(np.exp(-a))

def scaled_margin(X):
    z=X@B.T
    lr=log_rho(z)
    logs=np.column_stack([np.log(abs(COEFF[k]))+lr[:,i]+lr[:,j] for k,(i,j) in enumerate(PAIRS)])
    m=np.max(logs,axis=1)
    scaled=np.exp(logs-m[:,None])*np.sign(COEFF)[None,:]
    numerator=np.sum(scaled,axis=1)
    denominator=np.sum(np.abs(scaled),axis=1)
    return numerator/denominator,m+np.log(numerator)

def main():
    rng=np.random.default_rng(20260713)
    blocks=[]
    for R,n in [(5,200000),(20,200000),(100,200000),(500,200000),(2000,200000)]:
        X=rng.uniform(-R,R,size=(n,2))
        ratio,logJ=scaled_margin(X)
        k=int(np.argmin(ratio))
        assert np.all(ratio>0)
        blocks.append({
          'radius':R,'samples':n,
          'minimum_normalized_margin':float(ratio[k]),
          'log_J_at_margin_minimizer':float(logJ[k]),
          'point':X[k].tolist(),
        })
    report={'status':'nonrigorous stress test passed','blocks':blocks,'warning':'Floating point sampling is corroboration only. Log scaling prevents tail underflow; the analytic AGE proof is load bearing.'}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2)+'\n', newline="\n")
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
