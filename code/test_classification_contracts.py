#!/usr/bin/env python3
"""Regress missing winding faces and incomplete canonical certificate inventories."""
from copy import deepcopy
from fractions import Fraction
import json
import unittest

import check_alternating_classification_certificate as replay
import exact_alternating_classification as core


class ClassificationContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate=json.loads(replay.REPORT.read_text(encoding='utf-8'))

    def test_current_inventory(self):
        replay.validate_pattern_inventory(self.certificate)

    def test_omitted_pattern_rejected(self):
        malformed=deepcopy(self.certificate)
        malformed['patterns'].pop()
        with self.assertRaises(AssertionError):
            replay.validate_pattern_inventory(malformed)

    def test_repeated_pattern_with_unchanged_totals_rejected(self):
        malformed=deepcopy(self.certificate)
        self.assertEqual(malformed['patterns'][1]['chamber_count'],
                         malformed['patterns'][3]['chamber_count'])
        malformed['patterns'][1]=deepcopy(malformed['patterns'][3])
        self.assertEqual(len(malformed['patterns']),48)
        self.assertEqual(sum(r['chamber_count'] for r in malformed['patterns']),2448)
        with self.assertRaises(AssertionError):
            replay.validate_pattern_inventory(malformed)

    def test_noncanonical_same_chirotope_rejected(self):
        malformed=deepcopy(self.certificate)
        rec=malformed['patterns'][0]
        rec['directions'][0]=[2*x for x in rec['directions'][0]]
        D=rec['directions']
        signs=[1 if core.det2(D[i],D[j])>0 else -1 for i,j in core.PAIRS]
        self.assertEqual(signs,rec['pattern'])
        with self.assertRaises(AssertionError):
            replay.validate_pattern_inventory(malformed)

    def test_exact_missing_winding_three(self):
        directions=[(1,0),(-1,1),(0,-1),(1,1)]
        magnitudes=tuple(Fraction(x,28000000)
                         for x in (10769230,2692308,6461539,8076923))
        self.assertTrue(all(core.evalf(w,magnitudes)!=0
                            for w in core.walls_for(directions)))
        self.assertEqual(core.spectrum(directions,magnitudes),(0,1,2,3))
        edges=[(magnitudes[i]*d[0],magnitudes[i]*d[1])
               for i,d in enumerate(directions)]
        edges += [(-x,-y) for x,y in edges]
        vertices=[]
        current=(Fraction(0),Fraction(0))
        for edge in edges:
            vertices.append(current)
            current=core.p_add(current,edge)
        target=(Fraction(32307689,112000000),Fraction(10769229,112000000))
        self.assertEqual(core.wind(vertices,target),3)

    def test_arbitrarily_thin_face_has_no_iteration_cutoff(self):
        directions=[(1,0),(-1,1),(0,-1),(1,1)]
        interior=(Fraction(5,14),Fraction(3,28),Fraction(1,4),Fraction(2,7))
        known=tuple(Fraction(x,28000000)
                    for x in (10769230,2692308,6461539,8076923))
        boundary=tuple(a+(b-a)/(1-Fraction(1,10**6))
                       for a,b in zip(interior,known))
        magnitudes=tuple(a+(b-a)*(1-Fraction(1,10**100))
                         for a,b in zip(interior,boundary))
        self.assertTrue(all(core.evalf(w,magnitudes)!=0
                            for w in core.walls_for(directions)))
        self.assertEqual(core.spectrum(directions,magnitudes),(0,1,2,3))


if __name__=='__main__':
    unittest.main(verbosity=2)
