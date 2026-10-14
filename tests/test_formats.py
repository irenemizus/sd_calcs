import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import formats
import states

AU_TO_CM = 219474.624


def _hartree(cm):
    return cm / AU_TO_CM


class TestCalcFormat(unittest.TestCase):
    def _write_fort14(self, d, content):
        p = os.path.join(d, 'fort.14.test')
        with open(p, 'w') as f:
            f.write(content)
        return p

    def test_j0_two_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            content = (
                "0  0  0  0  0  2\n"
                "  {:.9f}  {:.9f}\n"
                "0  0  0  0  0  1\n"
                "  {:.9f}\n"
            ).format(_hartree(100.0), _hartree(200.0), _hartree(300.0))
            p = self._write_fort14(d, content)
            fmt = formats.CalcFormat(p, Ezero=0.0)
            result = fmt.parse_file()
            self.assertEqual(len(result), 3)
            self.assertEqual(result[0].J, 0)
            self.assertEqual(result[0].sym, 1)
            self.assertAlmostEqual(result[0].E, 100.0, delta=0.1)
            self.assertEqual(result[1].sym, 1)
            self.assertAlmostEqual(result[1].E, 200.0, delta=0.1)
            self.assertEqual(result[2].sym, 4)
            self.assertAlmostEqual(result[2].E, 300.0, delta=0.1)
            self.assertEqual(result[0].N, 1)
            self.assertEqual(result[1].N, 2)
            self.assertEqual(result[2].N, 1)

    def test_j1_four_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            content = (
                "1  0  0  0  0  1\n"
                "  {:.9f}\n"
                "1  0  0  0  0  1\n"
                "  {:.9f}\n"
                "1  0  0  0  0  1\n"
                "  {:.9f}\n"
                "1  0  0  0  0  1\n"
                "  {:.9f}\n"
            ).format(_hartree(100.0), _hartree(200.0), _hartree(300.0), _hartree(400.0))
            p = self._write_fort14(d, content)
            fmt = formats.CalcFormat(p, Ezero=0.0)
            result = fmt.parse_file()
            self.assertEqual(len(result), 4)
            self.assertEqual([s.sym for s in result], [1, 2, 3, 4])
            self.assertAlmostEqual(result[0].E, 100.0, delta=0.1)
            self.assertAlmostEqual(result[3].E, 400.0, delta=0.1)

    def test_e_zero_offset(self):
        with tempfile.TemporaryDirectory() as d:
            content = (
                "0  0  0  0  0  1\n"
                "  {:.9f}\n"
            ).format(_hartree(100.0))
            p = self._write_fort14(d, content)
            fmt = formats.CalcFormat(p, Ezero=50.0)
            result = fmt.parse_file()
            self.assertEqual(len(result), 1)
            self.assertAlmostEqual(result[0].E, 50.0, delta=0.1)

    def test_d_exponent(self):
        with tempfile.TemporaryDirectory() as d:
            content = (
                "0  0  0  0  0  1\n"
                "  {:.9f}D0\n"
            ).format(_hartree(100.0))
            p = self._write_fort14(d, content)
            fmt = formats.CalcFormat(p, Ezero=0.0)
            result = fmt.parse_file()
            self.assertEqual(len(result), 1)
            self.assertAlmostEqual(result[0].E, 100.0, delta=0.1)

    def test_multiple_energies_per_line(self):
        with tempfile.TemporaryDirectory() as d:
            content = (
                "0  0  0  0  0  4\n"
                "  {:.9f}  {:.9f}  {:.9f}  {:.9f}\n"
            ).format(_hartree(100.0), _hartree(200.0), _hartree(300.0), _hartree(400.0))
            p = self._write_fort14(d, content)
            fmt = formats.CalcFormat(p, Ezero=0.0)
            result = fmt.parse_file()
            self.assertEqual(len(result), 4)
            self.assertAlmostEqual(result[0].E, 100.0, delta=0.1)
            self.assertAlmostEqual(result[3].E, 400.0, delta=0.1)


class TestExpFormatH216O(unittest.TestCase):
    def _write_exp(self, d, lines):
        p = os.path.join(d, 'exp.txt')
        with open(p, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        return p

    def test_basic_parse(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "0 0 0 0 0 0 100.0",
                "0 0 0 0 0 0 200.0",
                "0 0 1 0 0 1 300.0",
            ]
            p = self._write_exp(d, lines)
            fmt = formats.ExpFormatH216O(p, J_list=['0'], J_place=3)
            result = fmt.parse_file()
            self.assertEqual(len(result), 3)
            self.assertEqual(result[0].J, 0)
            self.assertEqual(result[0].sym, 1)
            self.assertAlmostEqual(result[0].E, 100.0)
            self.assertEqual(result[1].sym, 1)
            self.assertAlmostEqual(result[1].E, 200.0)
            self.assertEqual(result[2].sym, 4)
            self.assertAlmostEqual(result[2].E, 300.0)

    def test_j_filtering(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "0 0 0 0 0 0 100.0",
                "0 0 1 1 0 1 300.0",
                "0 0 1 2 1 1 400.0",
            ]
            p = self._write_exp(d, lines)
            fmt = formats.ExpFormatH216O(p, J_list=['0', '1'], J_place=3)
            result = fmt.parse_file()
            self.assertEqual(len(result), 2)
            self.assertEqual([s.J for s in result], [0, 1])

    def test_sym_definition_j0(self):
        fmt = formats.ExpFormatH216O.__new__(formats.ExpFormatH216O)
        qn_a1 = states.QuantNumbersH216O(0, 0, 0, 0, 0, 0)
        qn_b2 = states.QuantNumbersH216O(0, 0, 1, 0, 0, 1)
        self.assertEqual(fmt.sym_definition(qn_a1), 1)
        self.assertEqual(fmt.sym_definition(qn_b2), 4)

    def test_sym_definition_j1(self):
        fmt = formats.ExpFormatH216O.__new__(formats.ExpFormatH216O)
        qn_b1 = states.QuantNumbersH216O(0, 0, 1, 1, 0, 1)
        self.assertEqual(fmt.sym_definition(qn_b1), 3)
        qn_a1 = states.QuantNumbersH216O(0, 0, 1, 1, 1, 1)
        self.assertEqual(fmt.sym_definition(qn_a1), 1)

    def test_qn_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            lines = ["0 0 1 1 0 1 300.0"]
            p = self._write_exp(d, lines)
            fmt = formats.ExpFormatH216O(p, J_list=['1'], J_place=3)
            result = fmt.parse_file()
            self.assertEqual(len(result), 1)
            qn = result[0].qn
            self.assertEqual(qn.v1, 0)
            self.assertEqual(qn.v2, 0)
            self.assertEqual(qn.v3, 1)
            self.assertEqual(qn.J, 1)
            self.assertEqual(qn.Ka, 0)
            self.assertEqual(qn.Kc, 1)


class TestExoMolFormatH216O(unittest.TestCase):
    def _write_states(self, d, lines):
        p = os.path.join(d, 'test.states')
        with open(p, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        return p

    def test_basic_parse(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "1  100.0  1  0  1",
                "2  200.0  1  0  1",
                "3  300.0  1  0  4",
            ]
            p = self._write_states(d, lines)
            fmt = formats.ExoMolFormatH216O(p, J_list=['0'], J_place=3)
            result = fmt.parse_file()
            self.assertEqual(len(result), 3)
            self.assertEqual(result[0].J, 0)
            self.assertEqual(result[0].sym, 1)
            self.assertAlmostEqual(result[0].E, 100.0)
            self.assertEqual(result[1].sym, 1)
            self.assertAlmostEqual(result[1].E, 200.0)
            self.assertEqual(result[2].sym, 4)
            self.assertAlmostEqual(result[2].E, 300.0)

    def test_N_values(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "1  100.0  1  0  1",
                "5  200.0  1  0  1",
            ]
            p = self._write_states(d, lines)
            fmt = formats.ExoMolFormatH216O(p, J_list=['0'], J_place=3)
            result = fmt.parse_file()
            self.assertEqual(result[0].N, 1)
            self.assertEqual(result[1].N, 5)

    def test_j_filtering(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "1  100.0  1  0  1",
                "2  200.0  1  1  1",
                "3  300.0  1  2  1",
            ]
            p = self._write_states(d, lines)
            fmt = formats.ExoMolFormatH216O(p, J_list=['0', '1'], J_place=3)
            result = fmt.parse_file()
            self.assertEqual(len(result), 2)
            self.assertEqual([s.J for s in result], [0, 1])

    def test_sym_definition_j0(self):
        fmt = formats.ExoMolFormatH216O.__new__(formats.ExoMolFormatH216O)
        qn = states.QuantNumbersJOnly(0)
        self.assertEqual(fmt.sym_definition(qn, lbl='1'), 1)
        self.assertEqual(fmt.sym_definition(qn, lbl='4'), 4)

    def test_sym_definition_j_odd(self):
        fmt = formats.ExoMolFormatH216O.__new__(formats.ExoMolFormatH216O)
        qn = states.QuantNumbersJOnly(1)
        self.assertEqual(fmt.sym_definition(qn, lbl='1'), 2)
        self.assertEqual(fmt.sym_definition(qn, lbl='2'), 3)
        self.assertEqual(fmt.sym_definition(qn, lbl='3'), 1)
        self.assertEqual(fmt.sym_definition(qn, lbl='4'), 4)

    def test_sym_definition_j_even(self):
        fmt = formats.ExoMolFormatH216O.__new__(formats.ExoMolFormatH216O)
        qn = states.QuantNumbersJOnly(2)
        self.assertEqual(fmt.sym_definition(qn, lbl='1'), 1)
        self.assertEqual(fmt.sym_definition(qn, lbl='2'), 4)
        self.assertEqual(fmt.sym_definition(qn, lbl='3'), 2)
        self.assertEqual(fmt.sym_definition(qn, lbl='4'), 3)


class TestExpFormatN2O(unittest.TestCase):
    def _write_exp(self, d, lines):
        p = os.path.join(d, 'exp_n2o.txt')
        with open(p, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        return p

    def test_basic_parse(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "0 1 0 e 0 100.0",
                "0 0 1 e 0 200.0",
            ]
            p = self._write_exp(d, lines)
            fmt = formats.ExpFormatN2O(p, J_list=['0'], J_place=4)
            result = fmt.parse_file()
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0].J, 0)
            self.assertEqual(result[0].sym, 1)
            self.assertAlmostEqual(result[0].E, 100.0)
            self.assertEqual(result[1].sym, 1)
            self.assertAlmostEqual(result[1].E, 200.0)

    def test_unc_and_nlines(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "0 1 0 e 0 100.0 0.001 10",
                "0 0 1 e 0 200.0 0.002 5",
            ]
            p = self._write_exp(d, lines)
            fmt = formats.ExpFormatN2O(p, J_list=['0'], J_place=4)
            result = fmt.parse_file()
            self.assertAlmostEqual(result[0].unc, 0.001)
            self.assertEqual(result[0].nlines, 10)
            self.assertAlmostEqual(result[1].unc, 0.002)
            self.assertEqual(result[1].nlines, 5)

    def test_no_unc_columns(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "0 1 0 e 0 100.0",
            ]
            p = self._write_exp(d, lines)
            fmt = formats.ExpFormatN2O(p, J_list=['0'], J_place=4)
            result = fmt.parse_file()
            self.assertEqual(result[0].unc, 0.0)
            self.assertEqual(result[0].nlines, 0)

    def test_sym_definition(self):
        fmt = formats.ExpFormatN2O.__new__(formats.ExpFormatN2O)
        self.assertEqual(fmt.sym_definition(lbl='e'), 1)
        self.assertEqual(fmt.sym_definition(lbl='f'), 2)

    def test_qn_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            lines = ["1 2 3 e 5 100.0"]
            p = self._write_exp(d, lines)
            fmt = formats.ExpFormatN2O(p, J_list=['5'], J_place=4)
            result = fmt.parse_file()
            qn = result[0].qn
            self.assertEqual(qn.v1, 1)
            self.assertEqual(qn.v2, 3)
            self.assertEqual(qn.l, 2)
            self.assertEqual(qn.J, 5)


class TestRITZFormatN2O556(unittest.TestCase):
    HEADER = "J isym Energy/cm-1   v1  v2  l2 v3 ief e/f unc/cm-1 #lines"

    def _write_ritz(self, d, lines):
        p = os.path.join(d, 'ritz.txt')
        all_lines = [self.HEADER] + lines
        with open(p, 'w') as f:
            f.write('\n'.join(all_lines) + '\n')
        return p

    def test_basic_parse(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "0 2  100.0  0  0  0  0  2  e  0.001  10",
                "0 1  200.0  0  1  0  0  2  e  0.001  10",
            ]
            p = self._write_ritz(d, lines)
            fmt = formats.RITZFormatN2O556(p, J_list=['0'], J_place=0)
            result = fmt.parse_file()
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0].J, 0)
            self.assertEqual(result[0].sym, 1)
            self.assertAlmostEqual(result[0].E, 100.0)
            self.assertEqual(result[1].sym, 2)
            self.assertAlmostEqual(result[1].E, 200.0)

    def test_header_filtered(self):
        with tempfile.TemporaryDirectory() as d:
            lines = ["0 2  100.0  0  0  0  0  2  e  0.001  10"]
            p = self._write_ritz(d, lines)
            fmt = formats.RITZFormatN2O556(p, J_list=['0'], J_place=0)
            result = fmt.parse_file()
            self.assertEqual(len(result), 1)

    def test_unc_and_nlines(self):
        with tempfile.TemporaryDirectory() as d:
            lines = ["0 2  100.0  0  0  0  0  2  e  0.003  20"]
            p = self._write_ritz(d, lines)
            fmt = formats.RITZFormatN2O556(p, J_list=['0'], J_place=0)
            result = fmt.parse_file()
            self.assertAlmostEqual(result[0].unc, 0.003)
            self.assertEqual(result[0].nlines, 20)

    def test_sym_definition_j0(self):
        fmt = formats.RITZFormatN2O556.__new__(formats.RITZFormatN2O556)
        qn = states.QuantNumbersN2O(0, 0, 0, 0, 0)
        self.assertEqual(fmt.sym_definition(qn, lbl='2'), 1)
        self.assertEqual(fmt.sym_definition(qn, lbl='1'), 2)

    def test_sym_definition_j_odd(self):
        fmt = formats.RITZFormatN2O556.__new__(formats.RITZFormatN2O556)
        qn = states.QuantNumbersN2O(0, 0, 1, 0, 0)
        self.assertEqual(fmt.sym_definition(qn, lbl='2'), 2)
        self.assertEqual(fmt.sym_definition(qn, lbl='1'), 1)

    def test_qn_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            lines = ["3 2  500.0  1  2  1  3  2  e  0.001  5"]
            p = self._write_ritz(d, lines)
            fmt = formats.RITZFormatN2O556(p, J_list=['3'], J_place=0)
            result = fmt.parse_file()
            qn = result[0].qn
            self.assertEqual(qn.J, 3)
            self.assertEqual(qn.v1, 1)
            self.assertEqual(qn.v2, 2)
            self.assertEqual(qn.l, 1)
            self.assertEqual(qn.v3, 3)

    def test_j_filtering(self):
        with tempfile.TemporaryDirectory() as d:
            lines = [
                "0 2  100.0  0  0  0  0  2  e  0.001  10",
                "1 2  200.0  0  1  0  0  2  e  0.001  10",
                "2 2  300.0  0  0  1  0  2  e  0.001  10",
            ]
            p = self._write_ritz(d, lines)
            fmt = formats.RITZFormatN2O556(p, J_list=['0', '1'], J_place=0)
            result = fmt.parse_file()
            self.assertEqual(len(result), 2)
            self.assertEqual([s.J for s in result], [0, 1])


if __name__ == '__main__':
    unittest.main()
