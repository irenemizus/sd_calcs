import io
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sdcalcs.states as states
from sdcalcs.comparison import do_comparison, exp_number_det
from sdcalcs.sd import (StandDev, _write_one_sd, assign_precision_weights, calc_sd,
                        outliers_by_N_out_rel, sd_calculation, sd_calculation_rel,
                        write_sds, write_sds_to_file, write_sds_to_file_rel)
from sdcalcs.utils import out_file_name_generation


def _st(E, J, sym, N=100, qn=None, **kw):
    if qn is None:
        qn = states.QuantNumbersJOnly(J)
    return states.State(E, J, sym, N=N, qn=qn, **kw)


def _qn6(v1=0, v2=0, v3=0, J=0, Ka=0, Kc=0):
    return states.QuantNumbersH216O(v1, v2, v3, J, Ka, Kc)


def _cs(E_exp, E_calc, J, sym, N, w, E_diff, status, qn=None, **kw):
    if qn is None:
        qn = _qn6(J=J)
    return states.ComparedState(E_calc, E_exp, J, sym, N, w, E_diff, status, qn, **kw)


class TestExpNumberDet(unittest.TestCase):
    def test_simple_match_found(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(200.0, 0, 1, N=101), _st(300.0, 1, 1, N=102)])
        exp = _st(100.1, 0, 1, qn=_qn6())
        sc_fit, cs = exp_number_det(calc, exp, eps=0.2, obs_calc_comp=True)
        self.assertIsNotNone(sc_fit)
        self.assertAlmostEqual(sc_fit.E, 100.0)
        self.assertIs(cs.status, states.Status.FOUND)
        self.assertAlmostEqual(cs.E_diff, 0.1, places=6)
        self.assertEqual(cs.N, 100)
        self.assertEqual(cs.w, 1.0)

    def test_best_of_multiple_candidates(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(100.05, 0, 1, N=101), _st(200.0, 0, 1, N=102), _st(300.0, 1, 1, N=103)])
        exp = _st(100.03, 0, 1, qn=_qn6())
        sc_fit, cs = exp_number_det(calc, exp, eps=0.2, obs_calc_comp=True)
        self.assertAlmostEqual(sc_fit.E, 100.05)
        self.assertIs(cs.status, states.Status.FOUND)
        self.assertAlmostEqual(cs.E_diff, -0.02, places=6)

    def test_outlier(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(300.0, 1, 1, N=101)])
        exp = _st(101.0, 0, 1, qn=_qn6())
        sc_fit, cs = exp_number_det(calc, exp, eps=0.2, obs_calc_comp=True)
        self.assertIsNotNone(sc_fit)
        self.assertIs(cs.status, states.Status.OUTLIER)
        self.assertEqual(cs.w, 0.0)
        self.assertAlmostEqual(cs.E_diff, 1.0, places=6)

    def test_not_found_no_match(self):
        calc = states.States([_st(100.0, 1, 1, N=100)])
        exp = _st(200.0, 0, 1, qn=_qn6())
        sc_fit, cs = exp_number_det(calc, exp, eps=0.2, obs_calc_comp=True)
        self.assertIsNone(sc_fit)
        self.assertIs(cs.status, states.Status.NOT_FOUND)
        self.assertEqual(cs.w, 0.0)
        self.assertEqual(cs.N, 0)

    def test_not_found_wrong_sym(self):
        calc = states.States([_st(100.0, 0, 1, N=100)])
        exp = _st(100.0, 0, 2, qn=_qn6())
        sc_fit, cs = exp_number_det(calc, exp, eps=0.2, obs_calc_comp=True)
        self.assertIsNone(sc_fit)
        self.assertIs(cs.status, states.Status.NOT_FOUND)

    def test_eps_boundary(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(300.0, 1, 1, N=101)])
        exp = _st(100.15, 0, 1, qn=_qn6())
        sc_fit, cs = exp_number_det(calc, exp, eps=0.2, obs_calc_comp=True)
        self.assertIs(cs.status, states.Status.FOUND)
        self.assertAlmostEqual(cs.E_diff, 0.15, places=6)

    def test_single_candidate(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(300.0, 1, 1, N=101)])
        exp = _st(100.1, 0, 1, qn=_qn6())
        sc_fit, cs = exp_number_det(calc, exp, eps=0.2, obs_calc_comp=True)
        self.assertAlmostEqual(sc_fit.E, 100.0)
        self.assertIs(cs.status, states.Status.FOUND)

    def test_prec_weight_propagated(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(300.0, 1, 1, N=101)])
        exp = _st(100.1, 0, 1, qn=_qn6(), prec_weight=0.5, is_hp=False)
        sc_fit, cs = exp_number_det(calc, exp, eps=0.2, obs_calc_comp=True)
        self.assertEqual(cs.prec_weight, 0.5)
        self.assertFalse(cs.is_hp)
        self.assertEqual(cs.w, 0.5)

    def test_obs_calc_comp_false_uses_E_exp(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(300.0, 1, 1, N=101)])
        qn = _qn6()
        exp = states.ComparedState(99.5, 100.1, 0, 1, 0, 0.0, 0.0, states.Status.FOUND, qn)
        sc_fit, cs = exp_number_det(calc, exp, eps=0.2, obs_calc_comp=False)
        self.assertAlmostEqual(cs.E_diff, 100.1 - 100.0, places=6)


class TestDoComparison(unittest.TestCase):
    def test_basic(self):
        calc = states.States([
            _st(100.0, 0, 1, N=100), _st(200.0, 0, 1, N=101),
            _st(300.0, 0, 2, N=102), _st(400.0, 1, 1, N=103),
        ])
        exp = states.States([
            _st(100.1, 0, 1, qn=_qn6(J=0)), _st(200.1, 0, 1, qn=_qn6(J=0)),
        ])
        result = do_comparison(calc, exp, eps=0.2)
        self.assertEqual(len(result), 2)
        self.assertIs(result[0].status, states.Status.FOUND)
        self.assertIs(result[1].status, states.Status.FOUND)
        self.assertAlmostEqual(result[0].E_diff, 0.1, places=6)
        self.assertAlmostEqual(result[1].E_diff, 0.1, places=6)

    def test_calc_not_modified(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(200.0, 0, 1, N=101), _st(300.0, 1, 1, N=102)])
        exp = states.States([_st(100.1, 0, 1, qn=_qn6())])
        do_comparison(calc, exp, eps=0.2)
        self.assertEqual(len(calc), 3)

    def test_N_propagated_to_exp(self):
        calc = states.States([
            _st(100.0, 0, 1, N=100), _st(200.0, 0, 1, N=101),
            _st(300.0, 0, 2, N=102), _st(400.0, 1, 1, N=103),
        ])
        exp = states.States([_st(100.1, 0, 1, qn=_qn6())])
        do_comparison(calc, exp, eps=0.2)
        self.assertEqual(exp[0].N, 100)

    def test_outlier(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(300.0, 1, 1, N=101)])
        exp = states.States([_st(105.0, 0, 1, qn=_qn6())])
        result = do_comparison(calc, exp, eps=0.2)
        self.assertEqual(len(result), 1)
        self.assertIs(result[0].status, states.Status.OUTLIER)
        self.assertEqual(result[0].w, 0.0)

    def test_not_found(self):
        calc = states.States([_st(100.0, 1, 1, N=100)])
        exp = states.States([_st(200.0, 0, 1, qn=_qn6())])
        result = do_comparison(calc, exp, eps=0.2)
        self.assertEqual(len(result), 0)

    def test_multiple_J(self):
        calc = states.States([
            _st(100.0, 0, 1, N=100), _st(200.0, 0, 1, N=101),
            _st(300.0, 1, 1, N=102), _st(400.0, 1, 1, N=103),
            _st(500.0, 2, 1, N=104),
        ])
        exp = states.States([
            _st(100.1, 0, 1, qn=_qn6(J=0)),
            _st(300.1, 1, 1, qn=_qn6(J=1)),
        ])
        result = do_comparison(calc, exp, eps=0.2)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0].J, 0)
        self.assertEqual(result[1].J, 1)

    def test_w_propagated_to_exp(self):
        calc = states.States([_st(100.0, 0, 1, N=100), _st(300.0, 1, 1, N=101)])
        exp = states.States([_st(105.0, 0, 1, qn=_qn6())])
        do_comparison(calc, exp, eps=0.2)
        self.assertEqual(exp[0].w, 0.0)


class TestSdCalculation(unittest.TestCase):
    def test_empty(self):
        sd_l, sd_m, sd_h = sd_calculation([], 15000, 25000)
        self.assertEqual(sd_l.Nsd, 0)
        self.assertEqual(sd_l.Nout, 0)
        self.assertEqual(sd_l.sd, 0.0)

    def test_all_found_low(self):
        cses = [
            _cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND),
            _cs(200.0, 199.8, 0, 1, 2, 1.0, -0.2, states.Status.FOUND),
        ]
        sd_l, sd_m, sd_h = sd_calculation(cses, 15000, 25000)
        self.assertEqual(sd_l.Nsd, 2)
        expected = math.sqrt((0.1**2 + 0.2**2) / 2)
        self.assertAlmostEqual(sd_l.sd, expected, places=6)

    def test_band_splitting(self):
        cses = [
            _cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND),
            _cs(20000.0, 20000.5, 0, 1, 2, 1.0, 0.5, states.Status.FOUND),
            _cs(30000.0, 30000.3, 0, 1, 3, 0.0, 0.3, states.Status.OUTLIER),
            _cs(30000.0, 30000.7, 0, 1, 4, 0.0, 0.7, states.Status.OUTLIER),
        ]
        sd_l, sd_m, sd_h = sd_calculation(cses, 15000, 25000)
        self.assertEqual(sd_l.Nsd, 1)
        self.assertAlmostEqual(sd_l.sd, 0.1, places=6)
        self.assertEqual(sd_m.Nsd, 1)
        self.assertAlmostEqual(sd_m.sd, 0.5, places=6)
        self.assertEqual(sd_h.Nsd, 0)
        self.assertEqual(sd_h.Nout, 2)
        self.assertAlmostEqual(sd_h.sd_out, math.sqrt((0.09 + 0.49) / 2), places=6)

    def test_boundary_E_tr_l(self):
        cs = _cs(15000.0, 15000.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND)
        sd_l, sd_m, sd_h = sd_calculation([cs], 15000, 25000)
        self.assertEqual(sd_l.Nsd, 1)

    def test_boundary_E_tr_h(self):
        cs = _cs(25000.0, 25000.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND)
        sd_l, sd_m, sd_h = sd_calculation([cs], 15000, 25000)
        self.assertEqual(sd_h.Nsd, 1)

    def test_weighted_sd(self):
        cses = [
            _cs(100.0, 100.1, 0, 1, 1, 0.5, 0.1, states.Status.FOUND),
            _cs(200.0, 200.1, 0, 1, 2, 0.5, 0.1, states.Status.FOUND),
        ]
        sd_l, _, _ = sd_calculation(cses, 15000, 25000)
        self.assertAlmostEqual(sd_l.sd, math.sqrt((0.01*0.5 + 0.01*0.5) / 2), places=6)


class TestAssignPrecisionWeights(unittest.TestCase):
    def _make_states(self, items):
        result = []
        for i, (unc, nlines) in enumerate(items):
            qn = _qn6(J=0)
            result.append(_st(float(i + 100), 0, 1, N=i+1, qn=qn, unc=unc, nlines=nlines))
        return states.States(result)

    def test_basic_split(self):
        st = self._make_states([(0.001, 1), (0.01, 1)])
        assign_precision_weights(st, dE_th=0.001)
        self.assertTrue(st[0].is_hp)
        self.assertAlmostEqual(st[0].prec_weight, 1.0)
        self.assertFalse(st[1].is_hp)
        self.assertAlmostEqual(st[1].prec_weight, 0.1, places=6)

    def test_hp_weight_param(self):
        st = self._make_states([(0.001, 1), (0.01, 1)])
        assign_precision_weights(st, dE_th=0.001, hp_weight=2.0)
        self.assertAlmostEqual(st[0].prec_weight, 2.0)

    def test_no_hp_levels(self):
        st = self._make_states([(0.01, 1), (0.02, 1)])
        assign_precision_weights(st, dE_th=0.001)
        self.assertFalse(st[0].is_hp)
        self.assertFalse(st[1].is_hp)
        self.assertAlmostEqual(st[0].prec_weight, 0.1, places=6)
        self.assertAlmostEqual(st[1].prec_weight, 0.05, places=6)

    def test_all_hp(self):
        st = self._make_states([(0.0005, 1), (0.001, 1)])
        assign_precision_weights(st, dE_th=0.001)
        self.assertTrue(st[0].is_hp)
        self.assertTrue(st[1].is_hp)
        self.assertAlmostEqual(st[0].prec_weight, 1.0)
        self.assertAlmostEqual(st[1].prec_weight, 1.0)

    def test_dE_calculation(self):
        st = self._make_states([(0.001, 4)])
        assign_precision_weights(st, dE_th=0.01)
        self.assertAlmostEqual(st[0].dE, 0.002, places=6)

    def test_dE_zero_nlines(self):
        st = self._make_states([(0.001, 0)])
        assign_precision_weights(st, dE_th=0.001)
        self.assertAlmostEqual(st[0].dE, 0.0)
        self.assertTrue(st[0].is_hp)

    def test_returns_dE_max(self):
        st = self._make_states([(0.001, 1), (0.01, 1)])
        dE_max = assign_precision_weights(st, dE_th=0.001)
        self.assertAlmostEqual(dE_max, 0.001, places=6)


class TestSdCalculationRel(unittest.TestCase):
    def test_six_values(self):
        cses = [
            _cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, is_hp=True),
            _cs(20000.0, 20000.1, 0, 1, 2, 0.5, 0.1, states.Status.FOUND, is_hp=False),
        ]
        result = sd_calculation_rel(cses, 15000, 25000)
        self.assertEqual(len(result), 6)

    def test_hp_low(self):
        cses = [
            _cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, is_hp=True),
            _cs(200.0, 200.2, 0, 1, 2, 1.0, 0.2, states.Status.FOUND, is_hp=True),
        ]
        result = sd_calculation_rel(cses, 15000, 25000)
        hp_l = result[0]
        self.assertEqual(hp_l.Nsd, 2)
        expected = math.sqrt((0.01 + 0.04) / 2)
        self.assertAlmostEqual(hp_l.sd, expected, places=6)

    def test_lp_middle(self):
        cses = [
            _cs(20000.0, 20000.3, 0, 1, 1, 0.5, 0.3, states.Status.FOUND, is_hp=False),
        ]
        result = sd_calculation_rel(cses, 15000, 25000)
        lp_m = result[3]
        self.assertEqual(lp_m.Nsd, 1)
        self.assertAlmostEqual(lp_m.sd, 0.3, places=6)

    def test_outlier_hp_high(self):
        cses = [
            _cs(30000.0, 30000.5, 0, 1, 1, 0.0, 0.5, states.Status.OUTLIER, is_hp=True),
        ]
        result = sd_calculation_rel(cses, 15000, 25000)
        hp_h = result[4]
        self.assertEqual(hp_h.Nout, 1)
        self.assertAlmostEqual(hp_h.sd_out, 0.5, places=6)

    def test_empty_groups(self):
        cses = [
            _cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, is_hp=True),
        ]
        result = sd_calculation_rel(cses, 15000, 25000)
        self.assertEqual(result[1].Nsd, 0)
        self.assertEqual(result[2].Nsd, 0)
        self.assertEqual(result[3].Nsd, 0)
        self.assertEqual(result[4].Nsd, 0)
        self.assertEqual(result[5].Nsd, 0)


class TestCalcSd(unittest.TestCase):
    def test_delegates_to_standard(self):
        cses = [_cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND)]
        result = calc_sd(cses, 15000, 25000, use_rel=False)
        self.assertEqual(len(result), 3)

    def test_delegates_to_rel(self):
        cses = [_cs(100.0, 100.1, 0, 1, 1, 1.0, 0.1, states.Status.FOUND, is_hp=True)]
        result = calc_sd(cses, 15000, 25000, use_rel=True)
        self.assertEqual(len(result), 6)


class TestOutliersByNOutRel(unittest.TestCase):
    def _make_cses(self):
        return [
            _cs(100.0, 100.01, 0, 1, 1, 1.0, 0.01, states.Status.FOUND),
            _cs(200.0, 200.02, 0, 1, 2, 1.0, 0.02, states.Status.FOUND),
            _cs(300.0, 300.03, 0, 1, 3, 1.0, 0.03, states.Status.FOUND),
            _cs(400.0, 400.04, 0, 1, 4, 1.0, 0.04, states.Status.FOUND),
        ]

    def test_N_out_rel_zero(self):
        result, eps_max = outliers_by_N_out_rel(self._make_cses(), 0.0)
        self.assertEqual(len(result), 4)
        for cs in result:
            self.assertIs(cs.status, states.Status.FOUND)

    def test_N_out_rel_50(self):
        result, eps_max = outliers_by_N_out_rel(self._make_cses(), 50.0)
        n_found = sum(1 for cs in result if cs.status == states.Status.FOUND)
        n_out = sum(1 for cs in result if cs.status == states.Status.OUTLIER)
        self.assertEqual(n_found, 2)
        self.assertEqual(n_out, 2)
        for cs in result:
            if cs.status == states.Status.OUTLIER:
                self.assertEqual(cs.w, 0.0)

    def test_N_out_rel_25(self):
        result, eps_max = outliers_by_N_out_rel(self._make_cses(), 25.0)
        n_found = sum(1 for cs in result if cs.status == states.Status.FOUND)
        n_out = sum(1 for cs in result if cs.status == states.Status.OUTLIER)
        self.assertEqual(n_found, 3)
        self.assertEqual(n_out, 1)

    def test_sorted_output(self):
        result, _ = outliers_by_N_out_rel(self._make_cses(), 50.0)
        keys = [(cs.J, cs.sym, cs.E_exp) for cs in result]
        self.assertEqual(keys, sorted(keys))

    def test_prec_weight_preserved(self):
        cses = [
            _cs(100.0, 100.01, 0, 1, 1, 0.7, 0.01, states.Status.FOUND, prec_weight=0.7),
            _cs(200.0, 200.02, 0, 1, 2, 0.7, 0.02, states.Status.FOUND, prec_weight=0.7),
        ]
        result, _ = outliers_by_N_out_rel(cses, 50.0)
        for cs in result:
            if cs.status == states.Status.FOUND:
                self.assertAlmostEqual(cs.w, 0.7)


class TestOutFileNameGeneration(unittest.TestCase):
    def test_with_extension(self):
        result = out_file_name_generation('/out', 'comp.txt', '5', 'H2-16O')
        self.assertIn('H2-16O', result)
        self.assertIn('J5', result)
        self.assertTrue(result.endswith('.txt'))

    def test_without_extension(self):
        result = out_file_name_generation('/out', 'comp', '5', 'H2-16O')
        self.assertIn('H2-16O', result)
        self.assertIn('J5', result)

    def test_J_suffix_with_eps(self):
        result = out_file_name_generation('/out', 'comp.txt', 'all_eps=0.2', 'N2O')
        self.assertIn('N2O', result)
        self.assertIn('all_eps=0.2', result)


class TestWriteSds(unittest.TestCase):
    def test_write_one_sd(self):
        sd = StandDev(Nsd=10, Nout=2, sd=0.1, sd_out=0.5)
        f = io.StringIO()
        _write_one_sd(f, 'l', sd)
        content = f.getvalue()
        self.assertIn('sd_l', content)
        self.assertIn('0.1000', content)
        self.assertIn('N_l = 10', content)

    def test_write_one_sd_empty(self):
        sd = StandDev()
        f = io.StringIO()
        _write_one_sd(f, 'l', sd)
        self.assertEqual(f.getvalue(), '')

    def test_write_sds_to_file(self):
        sds = [
            StandDev(Nsd=10, Nout=0, sd=0.1),
            StandDev(Nsd=5, Nout=1, sd=0.2, sd_out=0.5),
            StandDev(Nsd=0, Nout=0),
        ]
        f = io.StringIO()
        write_sds_to_file(f, sds)
        content = f.getvalue()
        self.assertIn('sd_l', content)
        self.assertIn('sd_m', content)
        self.assertNotIn('sd_h', content)

    def test_write_sds_to_file_rel(self):
        sds = [StandDev(Nsd=1) for _ in range(6)]
        f = io.StringIO()
        write_sds_to_file_rel(f, sds)
        content = f.getvalue()
        for label in ['hp_l', 'lp_l', 'hp_m', 'lp_m', 'hp_h', 'lp_h']:
            self.assertIn(label, content)

    def test_write_sds_dispatch(self):
        sds = [StandDev(Nsd=1) for _ in range(3)]
        f = io.StringIO()
        write_sds(f, sds, use_rel=False)
        self.assertIn('sd_l', f.getvalue())

        f2 = io.StringIO()
        sds6 = [StandDev(Nsd=1) for _ in range(6)]
        write_sds(f2, sds6, use_rel=True)
        self.assertIn('hp_l', f2.getvalue())


class TestStandDev(unittest.TestCase):
    def test_defaults(self):
        sd = StandDev()
        self.assertEqual(sd.Nsd, 0)
        self.assertEqual(sd.Nout, 0)
        self.assertEqual(sd.sd, 0.0)
        self.assertEqual(sd.sd_out, 0.0)

    def test_custom(self):
        sd = StandDev(Nsd=5, Nout=2, sd=0.1, sd_out=0.3)
        self.assertEqual(sd.Nsd, 5)
        self.assertEqual(sd.Nout, 2)
        self.assertEqual(sd.sd, 0.1)
        self.assertEqual(sd.sd_out, 0.3)


if __name__ == '__main__':
    unittest.main()
