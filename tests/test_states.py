import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import states


class TestState(unittest.TestCase):
    def test_defaults(self):
        s = states.State(E=100.0, J=0, sym=1)
        self.assertEqual(s.E, 100.0)
        self.assertEqual(s.J, 0)
        self.assertEqual(s.sym, 1)
        self.assertEqual(s.N, 0)
        self.assertEqual(s.w, 1.0)
        self.assertIsNone(s.qn)
        self.assertEqual(s.unc, 0.0)
        self.assertEqual(s.nlines, 0)
        self.assertEqual(s.dE, 0.0)
        self.assertTrue(s.is_hp)
        self.assertEqual(s.prec_weight, 1.0)

    def test_custom(self):
        qn = states.QuantNumbersH216O(0, 0, 1, 1, 0, 1)
        s = states.State(E=200.5, J=1, sym=3, N=5, w=0.5, qn=qn,
                         unc=0.001, nlines=10, dE=0.003,
                         is_hp=False, prec_weight=0.3)
        self.assertEqual(s.E, 200.5)
        self.assertEqual(s.J, 1)
        self.assertEqual(s.sym, 3)
        self.assertEqual(s.N, 5)
        self.assertEqual(s.w, 0.5)
        self.assertIs(s.qn, qn)
        self.assertEqual(s.unc, 0.001)
        self.assertEqual(s.nlines, 10)
        self.assertEqual(s.dE, 0.003)
        self.assertFalse(s.is_hp)
        self.assertEqual(s.prec_weight, 0.3)


class TestQuantNumbers(unittest.TestCase):
    def test_base_len_iter_getitem(self):
        qn = states.QuantNumbers([10, 20, 30])
        self.assertEqual(len(qn), 3)
        self.assertEqual(list(qn), [10, 20, 30])
        self.assertEqual(qn[0], 10)
        self.assertEqual(qn[2], 30)

    def test_h216o(self):
        qn = states.QuantNumbersH216O(v1=1, v2=2, v3=3, J=4, Ka=5, Kc=6)
        self.assertEqual(len(qn), 6)
        self.assertEqual(qn.v1, 1)
        self.assertEqual(qn.v2, 2)
        self.assertEqual(qn.v3, 3)
        self.assertEqual(qn.J, 4)
        self.assertEqual(qn.Ka, 5)
        self.assertEqual(qn.Kc, 6)
        self.assertEqual(list(qn), [1, 2, 3, 4, 5, 6])

    def test_n2o_4qn(self):
        qn = states.QuantNumbersN2O(v1=1, v2=2, J=3, l=4)
        self.assertEqual(len(qn), 4)
        self.assertEqual(qn.v1, 1)
        self.assertEqual(qn.v2, 2)
        self.assertEqual(qn.J, 3)
        self.assertEqual(qn.l, 4)
        self.assertIsNone(qn.v3)
        self.assertEqual(list(qn), [1, 2, 3, 4])

    def test_n2o_5qn(self):
        qn = states.QuantNumbersN2O(v1=1, v2=2, J=3, l=4, v3=5)
        self.assertEqual(len(qn), 5)
        self.assertEqual(qn.v3, 5)
        self.assertEqual(list(qn), [1, 2, 3, 4, 5])

    def test_j_only(self):
        qn = states.QuantNumbersJOnly(J=7)
        self.assertEqual(len(qn), 1)
        self.assertEqual(qn.J, 7)
        self.assertEqual(list(qn), [7])


class TestSymType(unittest.TestCase):
    def test_from_int(self):
        self.assertIs(states.SymType.from_int(1), states.SymType.A1)
        self.assertIs(states.SymType.from_int(2), states.SymType.A2)
        self.assertIs(states.SymType.from_int(3), states.SymType.B1)
        self.assertIs(states.SymType.from_int(4), states.SymType.B2)
        self.assertIs(states.SymType.from_int(5), states.SymType.UKN)

    def test_from_str(self):
        self.assertIs(states.SymType.from_str('A1'), states.SymType.A1)
        self.assertIs(states.SymType.from_str('a2'), states.SymType.A2)
        self.assertIs(states.SymType.from_str('B1'), states.SymType.B1)
        self.assertIs(states.SymType.from_str('b2'), states.SymType.B2)

    def test_values(self):
        self.assertEqual(states.SymType.A1.value, 1)
        self.assertEqual(states.SymType.A2.value, 2)
        self.assertEqual(states.SymType.B1.value, 3)
        self.assertEqual(states.SymType.B2.value, 4)
        self.assertEqual(states.SymType.UKN.value, 5)


class TestStatus(unittest.TestCase):
    def test_from_int(self):
        self.assertIs(states.Status.from_int(0), states.Status.FOUND)
        self.assertIs(states.Status.from_int(1), states.Status.OUTLIER)
        self.assertIs(states.Status.from_int(2), states.Status.NOT_FOUND)

    def test_names(self):
        self.assertEqual(states.Status.FOUND.name, 'FOUND')
        self.assertEqual(states.Status.OUTLIER.name, 'OUTLIER')
        self.assertEqual(states.Status.NOT_FOUND.name, 'NOT_FOUND')

    def test_values(self):
        self.assertEqual(states.Status.FOUND.value, 0)
        self.assertEqual(states.Status.OUTLIER.value, 1)
        self.assertEqual(states.Status.NOT_FOUND.value, 2)


class TestStates(unittest.TestCase):
    def _st(self, E, J, sym, qn=None):
        if qn is None:
            qn = states.QuantNumbersJOnly(J)
        return states.State(E, J, sym, qn=qn)

    def test_sorting_by_J_sym_E(self):
        s1 = self._st(100.0, 1, 1)
        s2 = self._st(50.0, 0, 1)
        s3 = self._st(75.0, 0, 2)
        s4 = self._st(25.0, 0, 1)
        coll = states.States([s1, s3, s2, s4])
        self.assertEqual([s.E for s in coll], [25.0, 50.0, 75.0, 100.0])
        self.assertEqual([s.J for s in coll], [0, 0, 0, 1])

    def test_len_iter_getitem(self):
        s1 = self._st(100.0, 0, 1)
        s2 = self._st(200.0, 1, 1)
        coll = states.States([s1, s2])
        self.assertEqual(len(coll), 2)
        self.assertEqual(len(list(coll)), 2)
        self.assertEqual(coll[0].E, 100.0)
        self.assertEqual(coll[1].E, 200.0)

    def test_empty(self):
        coll = states.States([])
        self.assertEqual(len(coll), 0)
        self.assertEqual(list(coll), [])

    def test_remove_item(self):
        s1 = self._st(100.0, 0, 1)
        s2 = self._st(200.0, 1, 1)
        coll = states.States([s1, s2])
        coll.remove_item(s1)
        self.assertEqual(len(coll), 1)
        self.assertEqual(coll[0].E, 200.0)

    def test_write_1qn(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = states.QuantNumbersJOnly(0)
            s = states.State(100.0, 0, 1, N=1, qn=qn)
            states.States([s]).write_to_file(p)
            with open(p) as f:
                parts = f.readline().strip().split()
            self.assertEqual(int(parts[0]), 0)
            self.assertEqual(parts[1], 'A1')
            self.assertEqual(int(parts[2]), 1)
            self.assertAlmostEqual(float(parts[3]), 100.0)

    def test_write_6qn(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = states.QuantNumbersH216O(v1=0, v2=0, v3=1, J=1, Ka=0, Kc=1)
            s = states.State(3151.6, 1, 3, N=2, w=1.0, qn=qn)
            states.States([s]).write_to_file(p)
            with open(p) as f:
                parts = f.readline().strip().split()
            self.assertEqual(int(parts[0]), 1)
            self.assertEqual(parts[1], 'B1')
            self.assertEqual(int(parts[2]), 2)
            self.assertAlmostEqual(float(parts[3]), 3151.6)
            self.assertEqual(int(parts[4]), 0)
            self.assertEqual(int(parts[5]), 0)
            self.assertEqual(int(parts[6]), 1)

    def test_write_4qn(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = states.QuantNumbersN2O(v1=0, v2=1, J=0, l=1)
            s = states.State(100.0, 0, 1, N=1, w=1.0, qn=qn)
            states.States([s]).write_to_file(p)
            with open(p) as f:
                parts = f.readline().strip().split()
            self.assertEqual(int(parts[0]), 0)
            self.assertEqual(parts[1], 'A1')
            self.assertAlmostEqual(float(parts[3]), 100.0)

    def test_write_5qn(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = states.QuantNumbersN2O(v1=0, v2=1, J=0, l=1, v3=2)
            s = states.State(100.0, 0, 1, N=1, w=1.0, qn=qn)
            states.States([s]).write_to_file(p)
            with open(p) as f:
                parts = f.readline().strip().split()
            self.assertEqual(int(parts[0]), 0)
            self.assertEqual(parts[1], 'A1')
            self.assertAlmostEqual(float(parts[3]), 100.0)


class TestComparedState(unittest.TestCase):
    def test_construction(self):
        qn = states.QuantNumbersH216O(0, 0, 1, 1, 0, 1)
        cs = states.ComparedState(3151.5, 3151.6, 1, 3, 2, 1.0, 0.1,
                                  states.Status.FOUND, qn)
        self.assertEqual(cs.E_calc, 3151.5)
        self.assertEqual(cs.E_exp, 3151.6)
        self.assertEqual(cs.E_diff, 0.1)
        self.assertEqual(cs.J, 1)
        self.assertEqual(cs.sym, 3)
        self.assertEqual(cs.N, 2)
        self.assertEqual(cs.w, 1.0)
        self.assertIs(cs.status, states.Status.FOUND)
        self.assertIs(cs.qn, qn)
        self.assertTrue(cs.is_hp)
        self.assertEqual(cs.prec_weight, 1.0)

    def test_defaults(self):
        cs = states.ComparedState(100.0, 100.1, 0, 1, 1, 1.0, 0.1)
        self.assertIsNone(cs.status)
        self.assertIsNone(cs.qn)
        self.assertTrue(cs.is_hp)
        self.assertEqual(cs.prec_weight, 1.0)


class TestComparisonList(unittest.TestCase):
    def _cs(self, E_exp, E_calc, J, sym, N, w, E_diff, status, qn):
        return states.ComparedState(E_calc, E_exp, J, sym, N, w, E_diff, status, qn)

    def _qn6(self, v1=0, v2=0, v3=0, J=0, Ka=0, Kc=0):
        return states.QuantNumbersH216O(v1, v2, v3, J, Ka, Kc)

    def _qn4(self, v1=0, v2=0, J=0, l=0):
        return states.QuantNumbersN2O(v1, v2, J, l)

    def test_sorting(self):
        c1 = self._cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, self._qn6())
        c2 = self._cs(200.0, 200.1, 0, 1, 2, 1.0, 0.1, states.Status.FOUND, self._qn6())
        c3 = self._cs(50.0, 50.1, 0, 2, 1, 1.0, 0.1, states.Status.FOUND, self._qn6())
        coll = states.ComparisonList([c3, c2, c1])
        self.assertEqual(coll[0].E_exp, 100.0)
        self.assertEqual(coll[1].E_exp, 200.0)
        self.assertEqual(coll[2].E_exp, 50.0)

    def test_add_item_reinserts_sorted(self):
        c1 = self._cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, self._qn6())
        c2 = self._cs(50.0, 50.1, 0, 1, 2, 1.0, 0.1, states.Status.FOUND, self._qn6())
        coll = states.ComparisonList([c1])
        coll.add_item(c2)
        self.assertEqual(coll[0].E_exp, 50.0)
        self.assertEqual(coll[1].E_exp, 100.0)

    def test_len_iter_getitem(self):
        c1 = self._cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, self._qn6())
        c2 = self._cs(200.0, 200.1, 0, 1, 2, 1.0, 0.1, states.Status.FOUND, self._qn6())
        coll = states.ComparisonList([c1, c2])
        self.assertEqual(len(coll), 2)
        self.assertEqual(len(list(coll)), 2)
        self.assertEqual(coll[0].E_exp, 100.0)

    def test_empty(self):
        coll = states.ComparisonList([])
        self.assertEqual(len(coll), 0)

    def test_write_6qn(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = self._qn6(v1=0, v2=0, v3=1, J=1, Ka=0, Kc=1)
            cs = self._cs(3151.6, 3151.5, 1, 3, 2, 1.0, 0.1, states.Status.FOUND, qn)
            states.ComparisonList([cs]).write_to_file(p)
            with open(p) as f:
                parts = f.readline().strip().split()
            self.assertEqual(len(parts), 14)
            self.assertEqual(int(parts[0]), 1)
            self.assertEqual(parts[1], 'B1')
            self.assertAlmostEqual(float(parts[3]), 3151.6)
            self.assertAlmostEqual(float(parts[4]), 3151.5)
            self.assertAlmostEqual(float(parts[5]), 0.1)
            self.assertEqual(parts[13], 'FOUND')

    def test_write_4qn(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = self._qn4(v1=0, v2=1, J=0, l=1)
            cs = self._cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, qn)
            states.ComparisonList([cs]).write_to_file(p)
            with open(p) as f:
                parts = f.readline().strip().split()
            self.assertEqual(len(parts), 12)
            self.assertEqual(int(parts[0]), 0)
            self.assertEqual(parts[1], 'A1')
            self.assertEqual(parts[11], 'FOUND')

    def test_write_filter_by_zero(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = self._qn6()
            c_found = self._cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, qn)
            c_out = self._cs(200.0, 200.5, 0, 1, 2, 0.0, 0.5, states.Status.OUTLIER, qn)
            result = states.ComparisonList([c_found, c_out]).write_to_file(p, filter_by=0.0)
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0].E_exp, 200.0)
            with open(p) as f:
                lines = f.readlines()
            self.assertEqual(len(lines), 1)

    def test_write_no_filter(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = self._qn6()
            c1 = self._cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, qn)
            c2 = self._cs(200.0, 200.5, 0, 1, 2, 0.0, 0.5, states.Status.OUTLIER, qn)
            result = states.ComparisonList([c1, c2]).write_to_file(p)
            self.assertEqual(len(result), 2)

    def test_parse_6qn(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = self._qn6(v1=0, v2=0, v3=1, J=1, Ka=0, Kc=1)
            cs = self._cs(3151.6, 3151.5, 1, 3, 2, 1.0, 0.1, states.Status.FOUND, qn)
            states.ComparisonList([cs]).write_to_file(p)
            parsed = states.ComparisonList([]).parse_file(p, 'exp.h2-16o.fmt')
            self.assertEqual(len(parsed), 1)
            self.assertAlmostEqual(parsed[0].E_exp, 3151.6)
            self.assertAlmostEqual(parsed[0].E_calc, 3151.5)
            self.assertAlmostEqual(parsed[0].E_diff, 0.1)
            self.assertEqual(parsed[0].J, 1)
            self.assertEqual(parsed[0].N, 2)
            self.assertAlmostEqual(parsed[0].w, 1.0)
            self.assertIs(parsed[0].status, states.Status.FOUND)
            self.assertEqual(parsed[0].qn.v3, 1)

    def test_parse_4qn(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = self._qn4(v1=0, v2=1, J=0, l=1)
            cs = self._cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, qn)
            states.ComparisonList([cs]).write_to_file(p)
            parsed = states.ComparisonList([]).parse_file(p, 'exp.n2o.fmt')
            self.assertEqual(len(parsed), 1)
            self.assertAlmostEqual(parsed[0].E_exp, 100.0)
            self.assertEqual(parsed[0].qn.v1, 0)
            self.assertEqual(parsed[0].qn.v2, 1)

    def test_round_trip_6qn(self):
        with tempfile.TemporaryDirectory() as d:
            p1 = os.path.join(d, '1.txt')
            p2 = os.path.join(d, '2.txt')
            qns = [
                self._qn6(0, 0, 0, 0, 0, 0),
                self._qn6(0, 0, 1, 1, 0, 1),
                self._qn6(0, 1, 0, 2, 1, 1),
            ]
            cses = [
                self._cs(100.0, 100.05, 0, 1, 1, 1.0, 0.05, states.Status.FOUND, qns[0]),
                self._cs(200.0, 200.3, 1, 3, 1, 0.0, 0.3, states.Status.OUTLIER, qns[1]),
                self._cs(300.0, 299.9, 2, 4, 1, 1.0, -0.1, states.Status.FOUND, qns[2]),
            ]
            states.ComparisonList(cses).write_to_file(p1)
            parsed = states.ComparisonList([]).parse_file(p1, 'exp.h2-16o.fmt')
            states.ComparisonList(parsed).write_to_file(p2)
            with open(p1) as f1, open(p2) as f2:
                self.assertEqual(f1.read(), f2.read())

    def test_round_trip_4qn(self):
        with tempfile.TemporaryDirectory() as d:
            p1 = os.path.join(d, '1.txt')
            p2 = os.path.join(d, '2.txt')
            qns = [
                self._qn4(0, 1, 0, 1),
                self._qn4(1, 0, 1, 0),
            ]
            cses = [
                self._cs(100.0, 100.05, 0, 1, 1, 1.0, 0.05, states.Status.FOUND, qns[0]),
                self._cs(200.0, 200.3, 1, 2, 1, 0.0, 0.3, states.Status.OUTLIER, qns[1]),
            ]
            states.ComparisonList(cses).write_to_file(p1)
            parsed = states.ComparisonList([]).parse_file(p1, 'exp.n2o.fmt')
            states.ComparisonList(parsed).write_to_file(p2)
            with open(p1) as f1, open(p2) as f2:
                self.assertEqual(f1.read(), f2.read())


class TestHITRANStates(unittest.TestCase):
    def test_sorting_by_J_sym_E0(self):
        qn1 = states.QuantNumbersH216O(0, 0, 0, 0, 0, 0)
        qn2 = states.QuantNumbersH216O(0, 0, 1, 0, 0, 1)
        s1 = states.State([200.0], 0, 1, qn=qn1)
        s2 = states.State([100.0], 0, 1, qn=qn2)
        coll = states.HITRANStates([s1, s2])
        self.assertEqual(coll[0].E[0], 100.0)
        self.assertEqual(coll[1].E[0], 200.0)

    def test_write_with_eth(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn1 = states.QuantNumbersH216O(0, 0, 0, 0, 0, 0)
            qn2 = states.QuantNumbersH216O(0, 0, 1, 0, 0, 1)
            s1 = states.State([100.0], 0, 1, qn=qn1)
            s2 = states.State([60000.0], 0, 1, qn=qn2)
            states.HITRANStates([s1, s2]).write_to_file(p, E_th=50000.0)
            with open(p) as f:
                content = f.read()
            self.assertIn('100.000000', content)
            self.assertNotIn('60000.000000', content)

    def test_write_multi_energy(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'out.txt')
            qn = states.QuantNumbersH216O(0, 0, 0, 0, 0, 0)
            s = states.State([100.0, 200.0], 0, 1, qn=qn)
            states.HITRANStates([s]).write_to_file(p)
            with open(p) as f:
                lines = f.readlines()
            self.assertEqual(len(lines), 2)


if __name__ == '__main__':
    unittest.main()
