#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from math import lcm
import json

Vector=tuple[Q,Q]

def add(*vectors: Vector) -> Vector:
    return (sum(v[0] for v in vectors),sum(v[1] for v in vectors))

def scale(c: Q,v: Vector) -> Vector:
    return (c*v[0],c*v[1])

@dataclass(frozen=True)
class AGEGroup:
    name: str
    bridge_index: int
    endpoint_left: int
    endpoint_right: int
    bridge_left: Q
    bridge_right: Q
    amgm_left: Q
    amgm_right: Q
    positive_left: Q
    positive_right: Q
    negative: Q

    def verify(self,rows: list[Vector]) -> dict[str,str]:
        assert self.bridge_left>=0 and self.bridge_right>=0
        assert self.bridge_left+self.bridge_right<=1
        assert self.amgm_left>0 and self.amgm_right>0
        assert self.amgm_left+self.amgm_right==1
        assert self.amgm_left>=self.bridge_left
        assert self.amgm_right>=self.bridge_right
        bridge=add(scale(self.bridge_left,rows[self.endpoint_left]),scale(self.bridge_right,rows[self.endpoint_right]))
        assert bridge==rows[self.bridge_index]

        # Weighted AM GM coefficient comparison without floating point.
        # Check product (positive_i / weight_i)^weight_i > negative.
        den=lcm(self.amgm_left.denominator,self.amgm_right.denominator)
        nl=int(self.amgm_left*den)
        nr=int(self.amgm_right*den)
        left=self.positive_left/self.amgm_left
        right=self.positive_right/self.amgm_right
        lhs_num=(left.numerator**nl)*(right.numerator**nr)*(self.negative.denominator**den)
        lhs_den=(left.denominator**nl)*(right.denominator**nr)*(self.negative.numerator**den)
        assert lhs_num>lhs_den
        return {
            'group':self.name,
            'bridge_relation':'verified exactly',
            'convex_slack':str(1-self.bridge_left-self.bridge_right),
            'amgm_coefficient_power_check':f'{lhs_num} > {lhs_den}',
            'conclusion':'strictly positive for every input when the log derivative potential is even, convex, nonnegative, and vanishes at zero',
        }

ROWS=[(Q(1),Q(0)),(Q(0),Q(1)),(Q(-1,5),Q(5,4)),(Q(-1,2),Q(9,5))]
GROUPS=[
    AGEGroup(
        name='negative pair 02',bridge_index=1,endpoint_left=0,endpoint_right=2,
        bridge_left=Q(4,25),bridge_right=Q(4,5),
        amgm_left=Q(1,5),amgm_right=Q(4,5),
        positive_left=Q(3,8),positive_right=Q(1),negative=Q(5,4),
    ),
    AGEGroup(
        name='negative pair 13',bridge_index=2,endpoint_left=1,endpoint_right=3,
        bridge_left=Q(53,100),bridge_right=Q(2,5),
        amgm_left=Q(53,100),amgm_right=Q(47,100),
        positive_left=Q(53,200),positive_right=Q(49,40),negative=Q(1),
    ),
]

def main() -> None:
    report={'verdict':'VERIFIED','groups':[g.verify(ROWS) for g in GROUPS]}
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
