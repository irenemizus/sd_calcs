import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN_PY = os.path.join(PROJ_ROOT, 'main.py')
PYTHON = sys.executable

AU_TO_CM = 219474.624


def _hartree(cm):
    return cm / AU_TO_CM


def _run_main(workdir, *args):
    cmd = [PYTHON, MAIN_PY] + list(args)
    result = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=30)
    return result


def _read_file(path):
    with open(path) as f:
        return f.read()


def _setup_h2o_exomol(workdir):
    inp = os.path.join(workdir, 'input', 'H2-16O')
    os.makedirs(inp, exist_ok=True)
    os.makedirs(os.path.join(workdir, 'output'), exist_ok=True)
    states_content = (
        "1  100.0  1  0  1\n"
        "2  200.0  1  0  1\n"
        "3  300.0  1  0  4\n"
    )
    with open(os.path.join(inp, 'test.states'), 'w') as f:
        f.write(states_content)
    exp_content = (
        "0 0 0 0 0 0 100.0\n"
        "0 0 0 0 0 0 200.0\n"
        "0 0 1 0 0 1 300.0\n"
    )
    with open(os.path.join(inp, 'exp.txt'), 'w') as f:
        f.write(exp_content)


def _setup_n2o_rel(workdir):
    inp = os.path.join(workdir, 'input', 'N2O')
    os.makedirs(inp, exist_ok=True)
    os.makedirs(os.path.join(workdir, 'output'), exist_ok=True)
    fort_content = (
        "0  0  0  0  0  4\n"
        "  {:.9f}  {:.9f}  {:.9f}  {:.9f}\n"
    ).format(_hartree(100.0), _hartree(200.0), _hartree(300.0), _hartree(400.0))
    with open(os.path.join(inp, 'fort.14.0'), 'w') as f:
        f.write(fort_content)
    # dE = unc * sqrt(nlines); HP if dE <= dE_th (0.001)
    # unc=0.0001, nlines=1 -> dE=0.0001 (HP); unc=0.01, nlines=10 -> dE=0.0316 (LP)
    exp_content = (
        "0 1 0 e 0 100.0 0.0001 1\n"
        "0 1 1 e 0 200.0 0.01 10\n"
        "0 1 2 e 0 300.0 0.0001 1\n"
        "0 1 3 e 0 400.0 0.01 10\n"
    )
    with open(os.path.join(inp, 'exp_n2o.txt'), 'w') as f:
        f.write(exp_content)


def _setup_n2o_basic(workdir):
    inp = os.path.join(workdir, 'input', 'N2O')
    os.makedirs(inp, exist_ok=True)
    os.makedirs(os.path.join(workdir, 'output'), exist_ok=True)
    fort_content = (
        "0  0  0  0  0  2\n"
        "  {:.9f}  {:.9f}\n"
    ).format(_hartree(100.0), _hartree(200.0))
    with open(os.path.join(inp, 'fort.14.0'), 'w') as f:
        f.write(fort_content)
    exp_content = (
        "0 1 0 e 0 100.0\n"
        "0 0 1 e 0 200.0\n"
    )
    with open(os.path.join(inp, 'exp_n2o.txt'), 'w') as f:
        f.write(exp_content)


def _setup_suboptimal(workdir):
    inp = os.path.join(workdir, 'input', 'H2-16O')
    os.makedirs(inp, exist_ok=True)
    os.makedirs(os.path.join(workdir, 'output'), exist_ok=True)
    states_content = (
        "100  100.0  1  0  1\n"
        "101  100.1  1  0  1\n"
        "102  200.0  1  0  1\n"
        "103  300.0  1  0  4\n"
        "104  400.0  1  1  1\n"
    )
    with open(os.path.join(inp, 'test.states'), 'w') as f:
        f.write(states_content)
    exp_content = (
        "0 0 0 0 0 0 100.05\n"
        "0 0 0 0 0 0 100.15\n"
        "0 0 0 0 0 0 200.0\n"
        "0 0 1 0 0 1 300.0\n"
    )
    with open(os.path.join(inp, 'exp.txt'), 'w') as f:
        f.write(exp_content)


class TestE2EH2OExoMol(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        _setup_h2o_exomol(self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_basic_run(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'test.states',
                      '--mol_name', 'H2-16O',
                      '--file_exp_name', 'exp.txt')
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        out_dir = os.path.join(self.tmpdir, 'output', 'H2-16O')
        comp_files = [f for f in os.listdir(out_dir) if f.startswith('comp')]
        self.assertTrue(len(comp_files) >= 1)
        comp_file = os.path.join(out_dir, comp_files[0])
        content = _read_file(comp_file)
        self.assertIn('100.0', content)
        self.assertIn('200.0', content)
        self.assertIn('300.0', content)

    def test_ens_file_generated(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'test.states',
                      '--mol_name', 'H2-16O',
                      '--file_exp_name', 'exp.txt')
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        out_dir = os.path.join(self.tmpdir, 'output', 'H2-16O')
        ens_files = [f for f in os.listdir(out_dir) if f.startswith('ens')]
        self.assertTrue(len(ens_files) >= 1)

    def test_sd_file_generated(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'test.states',
                      '--mol_name', 'H2-16O',
                      '--file_exp_name', 'exp.txt')
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        out_dir = os.path.join(self.tmpdir, 'output', 'H2-16O')
        sd_files = [f for f in os.listdir(out_dir) if 'sd' in f]
        self.assertTrue(len(sd_files) >= 1)


class TestE2EH2OExoMolECut(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        _setup_h2o_exomol(self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_E_cut(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'test.states',
                      '--mol_name', 'H2-16O',
                      '--file_exp_name', 'exp.txt',
                      '--E_cut', '250.0')
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        out_dir = os.path.join(self.tmpdir, 'output', 'H2-16O')
        ens_files = [f for f in os.listdir(out_dir) if f.startswith('ens')]
        self.assertTrue(len(ens_files) >= 1)
        ens_content = _read_file(os.path.join(out_dir, ens_files[0]))
        lines = [l for l in ens_content.strip().split('\n') if l.strip()]
        for line in lines:
            parts = line.split()
            if len(parts) >= 7:
                E = float(parts[6])
                self.assertLessEqual(E, 250.0, msg=f"Level with E={E} exceeds E_cut=250.0")


class TestE2EN2OBasic(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        _setup_n2o_basic(self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_basic_run(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'all',
                      '--mol_name', 'N2O',
                      '--file_exp_name', 'exp_n2o.txt')
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        out_dir = os.path.join(self.tmpdir, 'output', 'N2O')
        comp_files = [f for f in os.listdir(out_dir) if f.startswith('comp')]
        self.assertTrue(len(comp_files) >= 1)
        comp_content = _read_file(os.path.join(out_dir, comp_files[0]))
        self.assertIn('100.0', comp_content)


class TestE2EN2ORelWeights(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        _setup_n2o_rel(self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_rel_weights_run(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'all',
                      '--mol_name', 'N2O',
                      '--file_exp_name', 'exp_n2o.txt',
                      '--rel_weights',
                      '--dE_th', '0.001')
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        out_dir = os.path.join(self.tmpdir, 'output', 'N2O')
        sd_files = [f for f in os.listdir(out_dir) if 'sd' in f]
        self.assertTrue(len(sd_files) >= 1)
        found_rel = False
        for sf in sd_files:
            content = _read_file(os.path.join(out_dir, sf))
            if 'hp_l' in content:
                found_rel = True
                self.assertIn('lp_l', content)
        self.assertTrue(found_rel, "No sd file with rel-weight sections found")

    def test_hp_weight(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'all',
                      '--mol_name', 'N2O',
                      '--file_exp_name', 'exp_n2o.txt',
                      '--rel_weights',
                      '--dE_th', '0.001',
                      '--hp_weight', '2.0')
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        self.assertIn('high-precision', r.stdout)

    def test_rel_weights_requires_make_comp_true(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'all',
                      '--mol_name', 'N2O',
                      '--file_exp_name', 'exp_n2o.txt',
                      '--rel_weights',
                      '--make_comp_files', 'False')
        self.assertNotEqual(r.returncode, 0)


class TestE2ENOutRel(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        _setup_h2o_exomol(self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_N_out_rel(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'test.states',
                      '--mol_name', 'H2-16O',
                      '--file_exp_name', 'exp.txt',
                      '--N_out_rel', '33.0')
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        out_dir = os.path.join(self.tmpdir, 'output', 'H2-16O')
        comp_files = [f for f in os.listdir(out_dir) if f.startswith('comp')]
        self.assertTrue(len(comp_files) >= 1)


class TestE2EMakeCompFilesFalse(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        _setup_n2o_basic(self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_make_comp_false_after_true(self):
        r1 = _run_main(self.tmpdir,
                       '--mode', 'all',
                       '--mol_name', 'N2O',
                       '--file_exp_name', 'exp_n2o.txt',
                       '--make_comp_files', 'True')
        self.assertEqual(r1.returncode, 0, msg=r1.stderr)
        out_dir = os.path.join(self.tmpdir, 'output', 'N2O')
        files_after_run1 = os.listdir(out_dir)

        r2 = _run_main(self.tmpdir,
                       '--mode', 'all',
                       '--mol_name', 'N2O',
                       '--file_exp_name', 'exp_n2o.txt',
                       '--make_comp_files', 'False')
        self.assertEqual(r2.returncode, 0, msg=r2.stderr)
        files_after_run2 = os.listdir(out_dir)
        new_files = set(files_after_run2) - set(files_after_run1)
        self.assertTrue(any('sd' in f for f in new_files),
                        f"No new sd files after run2. files_after_run1={files_after_run1}, files_after_run2={files_after_run2}, stdout={r2.stdout}, stderr={r2.stderr}")


class TestE2ESuboptimal(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        _setup_suboptimal(self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_suboptimal_matches_file(self):
        r = _run_main(self.tmpdir,
                      '--mode', 'test.states',
                      '--mol_name', 'H2-16O',
                      '--file_exp_name', 'exp.txt')
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        out_dir = os.path.join(self.tmpdir, 'output', 'H2-16O')
        all_files = os.listdir(out_dir)
        self.assertIn('suboptimal_matches.txt', all_files,
                      f"suboptimal_matches.txt not found in {all_files}. stdout={r.stdout}")


if __name__ == '__main__':
    unittest.main()
