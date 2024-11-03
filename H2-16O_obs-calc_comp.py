import argparse
import copy
import math
import os

import formats
import main
import states


class ComparedState:
    def __init__(self, E_calc, E_exp, E_calc_old, sym, J, N, w, E_diff, E_diff_old, diff_rel, mark='', status=None, qn=None):
        self.E_calc = E_calc
        self.E = E_exp
        self.E_calc_old = E_calc_old
        self.E_diff = E_diff
        self.E_diff_old = E_diff_old
        self.diff_rel = diff_rel
        self.J = J
        self.sym = sym
        self.N = N
        self.w = w
        self.status = status
        self.mark = mark
        self.qn = qn


class ComparisonList:
    def __init__(self, list_comp_states):
        self.__comp_states = list_comp_states
        self.__comp_states = self.__sort()   # always sorted

    def __len__(self):
        return len(self.__comp_states)

    def __iter__(self):
        return iter(self.__comp_states)

    def __getitem__(self, item):
        return self.__comp_states[item]

    def add_item(self, item):
        self.__comp_states.append(item)
        self.__comp_states = self.__sort()

    def __sort(self):
        return sorted(self.__comp_states, key=lambda x: (x.J, x.sym, x.E))

    def write_to_file(self, out_file_name, filter_by=-1.0, mode='w'):
        with open(out_file_name, mode) as f:
            if filter_by < 0.0:
                states_to_write = self.__comp_states
            else:
                # write only the states with w = filter_by
                states_to_write = list(filter(lambda x: abs(x.w - filter_by) < 1e-3, self.__comp_states))

            for state in states_to_write:
                f.write("{0:2d} {1:3s} {2:5d} {3:16.6f} {4:16.6f} {5:16.6f} {6:4d} {7:3d} {8:3d} {9:3d} {10:3d} {11:3d} "
                        "{12:16.6f} {13:16.6f} {14:13.3f} {15:10s}\n".format(
                    state.J, states.SymType(state.sym).name, state.N, state.E, state.E_calc, state.E_diff,
                    state.qn.v1, state.qn.v2, state.qn.v3, state.qn.J, state.qn.Ka, state.qn.Kc, state.E_calc_old,
                    state.E_diff_old, state.diff_rel, state.mark
                ))

        return states_to_write


def do_comparison(states_comp, states_calc_old, comp_by_qns, eps=0.0, obs_calc_comp=False):
    states_old_cpy = copy.deepcopy(states_calc_old)
    comp_states = []

    # iterating by compared data
    for sc in states_comp:
        if comp_by_qns:
            for so in states_old_cpy:
                if sc.J == so.J and sc.qn.Ka == so.qn.Ka and sc.qn.Kc == so.qn.Kc and sc.qn.v1 == so.qn.v1 and \
                     sc.qn.v2 == so.qn.v2 and sc.qn.v3 == so.qn.v3:
                    E_diff_list = []
                    for E_old in so.E:
                        E_diff_list.append(sc.E_exp - E_old)
                    E_diff_abs_list = map(lambda x: abs(x), E_diff_list)
                    argmin = min(enumerate(E_diff_abs_list), key=lambda x: x[1])[0]   # argmin(E_diff_abs_list)
                    E_diff_old = E_diff_list[argmin]
                    E_old_best = so.E[argmin]
                    diff_rel = E_diff_old if sc.E_diff == 0.0 else E_diff_old / sc.E_diff
                    mark = "!!!" if abs(diff_rel) > 0.5 else ""
                    comp_old_state = ComparedState(sc.E_calc, sc.E_exp, E_old_best, sc.sym, sc.J, sc.N, sc.w,
                                                 sc.E_diff, E_diff_old, diff_rel, mark, sc.status, sc.qn)
                    comp_states.append(comp_old_state)
                    states_old_cpy.remove_item(so)
                    break
        else:
            # searching for the best fit among the old calculated energies
            sc_fit, compared_state = main.exp_number_det(states_old_cpy, sc, eps, obs_calc_comp)
            if sc_fit:
                states_old_cpy.remove_item(sc_fit)  # removing the found calculated energy
                E_old_best = compared_state.E_calc
                E_diff_old = compared_state.E_diff
                diff_rel = E_diff_old if sc.E_diff == 0.0 else E_diff_old / sc.E_diff
                mark = "!!!" if abs(diff_rel) > 0.5 else ""
                comp_old_state = ComparedState(sc.E_calc, sc.E_exp, E_old_best, sc.sym, sc.J, sc.N, sc.w,
                                               sc.E_diff, E_diff_old, diff_rel, mark, sc.status, sc.qn)
                comp_states.append(comp_old_state)

    return ComparisonList(comp_states), states_old_cpy


if __name__ == "__main__":
    """
        A script for finding related pairs of previous (HITRAN/POKAZATEL) and current calculated and observed states for H2-16O 
        molecule and obtaining deviations for them.

        Should be considered as a small separated add-on to the main sd_calcs script.
    """

    # Usually non-changeable parameters
    input_folder = "input"                  # A folder containing input data
    output_folder = "output"                # A folder containing output data

    # Command line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("--mol_name", type=str, help="name of the molecule", choices=['H2-16O'],
                        default='H2-16O')
    parser.add_argument("--file_comp_name", type=str, help="name of the input file with pre-generated "
                                                           "'current calculated vs observed' comparison data")
    parser.add_argument("--file_old_calc_name", type=str, help="name of the input file with previous (HITRAN/POKAZATEL) calculated data")

    parser.add_argument("--out_file_comp_name", type=str,
                        help="basic part of the output, containing the entire comparison results, file name; will be extended with the current parameter values",
                        default="comp.txt")
    parser.add_argument("--eps", type=float,
                        help="calculated energies with |obs-calc| <= eps will be considered as the 'FOUND' ones; the others will be marked as 'OUTLIER's", default=0.2)
    parser.add_argument("--Jmax", type=int,
                        help="maximum value of J for comparison", default=100)


    args = parser.parse_args()

    mol_name = args.mol_name
    file_comp_name = args.file_comp_name
    file_old_calc_name = args.file_old_calc_name
    out_file_comp_name = args.out_file_comp_name
    eps = args.eps
    Jmax = args.Jmax

    # Checking input directory
    full_inp_folder = os.path.join(os.getcwd(), input_folder, mol_name + "_old_calc_comp")
    if not os.path.isdir(full_inp_folder):
        print("Input folder doesn't exist")
        exit(1)

    # Creating/checking output directories
    full_out_folder = os.path.join(os.getcwd(), output_folder, mol_name + "_old_calc_comp")
    is_out_dir_exists = os.path.isdir(full_out_folder)
    if not is_out_dir_exists:
        os.makedirs(full_out_folder)

    out_file_name = os.path.join(full_out_folder, out_file_comp_name)

    J_list = [ str(J) for J in range(Jmax + 1) ]

    # Getting the full path for the input file with comparison data
    file_comp_name_full = os.path.join(full_inp_folder, file_comp_name)
    # Getting the full path for the input file with HITRAN data
    file_old_calc_name_full = os.path.join(full_inp_folder, file_old_calc_name)

    comp_list = states.ComparisonList([])
    comp_states = comp_list.parse_file(file_comp_name_full)

    # For HITRAN data -------------------
    h_format = formats.HITRANFormatH216O(file_old_calc_name_full, J_list, 3)
    hitran_states = h_format.parse_file()

    out_HITRAN_levs_file_name = os.path.join(full_out_folder, 'HITRAN_levs.txt')
    hitran_states.write_to_file(out_HITRAN_levs_file_name)

    comp_h_states, h_states_nf = do_comparison(comp_states, hitran_states, True)
    comp_h_states.write_to_file(out_file_name)

    out_nf_HITRAN_levs_file_name = os.path.join(full_out_folder, 'HITRAN_levs_not_found.txt')
    h_states_nf.write_to_file(out_nf_HITRAN_levs_file_name, 26267.8)
    # -----------------------------------

    # For POKAZATEL data ----------------
    # pokaz_format = formats.POKAZATELFormatH216O(file_old_calc_name_full, J_list, 3)
    # pokaz_states = pokaz_format.parse_file()
    #
    # comp_pokaz_states, pokaz_states_nf = do_comparison(comp_states, pokaz_states, False, eps)
    # comp_pokaz_states.write_to_file(out_file_name)
    #
    # out_nf_POKAZ_levs_file_name = os.path.join(full_out_folder, 'POKAZ_levs_not_found.txt')
    # pokaz_states_nf.write_to_file(out_nf_POKAZ_levs_file_name)
    # -----------------------------------



