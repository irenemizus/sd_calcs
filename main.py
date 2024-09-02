import copy
import math
import os.path
import shutil
import sys

import formats
import states


def exp_number_det(states_calc, state_exp, eps):
    sign = lambda x: math.copysign(1, x)
    status = states.Status.NOT_FOUND
    list_diff = []
    indx_next = 0

    for sc in states_calc:
        indx_next += 1
        if state_exp.J == sc.J and state_exp.sym == sc.sym:
            sc_next = sc
            diff_oc_first = state_exp.E - sc_next.E
            list_diff.append((sc_next, diff_oc_first))
            break

    for i in range(indx_next, len(states_calc)):
        sc_next = states_calc[i]
        if state_exp.J == sc_next.J and state_exp.sym == sc_next.sym:
            diff_oc_next = state_exp.E - sc_next.E
            list_diff.append((sc_next, diff_oc_next))

            if sign(diff_oc_first) != sign(diff_oc_next) or abs(diff_oc_first) <= abs(diff_oc_next):
                break

    if sc_next.N == len(states_calc):
        return None, states.ComparedState(0.0, state_exp.E, state_exp.J, state_exp.sym, 0, 0.0, 0.0, status, state_exp.qn)

    if list_diff:
        (sc_fit, diff_oc_fit) = list_diff[-1] if abs(list_diff[-1][1]) <= abs(list_diff[-2][1]) else list_diff[-2]
        if abs(diff_oc_fit) <= eps:
            return sc_fit, states.ComparedState(sc_fit.E, state_exp.E, state_exp.J, state_exp.sym, sc_fit.N, 1.0, diff_oc_fit, states.Status.FOUND, state_exp.qn)
        elif abs(diff_oc_fit) > eps:
            return sc_fit, states.ComparedState(sc_fit.E, state_exp.E, state_exp.J, state_exp.sym, sc_fit.N, 0.0, diff_oc_fit, states.Status.OUTLIER, state_exp.qn)


def do_comparison(states_calc, states_exp, eps):
    states_calc_cpy = copy.deepcopy(states_calc)
    comp_states = []
    for se in states_exp:
        sc_fit, compared_state = exp_number_det(states_calc_cpy, se, eps)
        if sc_fit:
            states_calc_cpy.remove_item(sc_fit)
            se.N = compared_state.N
            se.w = compared_state.w
            comp_states.append(compared_state)
    return states.ComparisonList(comp_states)


class StandDev:
    def __init__(self, Nsd=0, Nout=0, sd=0.0):
        self.sd = sd
        self.Nsd = Nsd
        self.Nout = Nout


def sd_calculation(comp_states, E_tr_l, E_tr_h):
    sd_l = StandDev()
    sd_m = StandDev()
    sd_h = StandDev()

    for cs in comp_states:
        if cs.E_exp <= E_tr_l:
            if cs.w > 0:
                sd_l.Nsd += 1
                sd_l.sd += cs.E_diff ** 2 * cs.w
            else:
                sd_l.Nout += 1
        elif cs.E_exp >= E_tr_h:
            if cs.w > 0:
                sd_h.Nsd += 1
                sd_h.sd += cs.E_diff**2 * cs.w
            else:
                sd_h.Nout += 1
        else:
            if cs.w > 0:
                sd_m.Nsd += 1
                sd_m.sd += cs.E_diff**2 * cs.w
            else:
                sd_m.Nout += 1

    if sd_h.Nsd > 0:
        sd_h.sd = math.sqrt(sd_h.sd / sd_h.Nsd)
    if sd_m.Nsd > 0:
        sd_m.sd = math.sqrt(sd_m.sd / sd_m.Nsd)
    if sd_l.Nsd > 0:
        sd_l.sd = math.sqrt(sd_l.sd / sd_l.Nsd)

    return sd_l, sd_m, sd_h

def write_sds_to_file(fo, sd_list):
    suffix = ["l", "m", "h"]
    s = 0
    for sd in sd_list:
        if sd.Nsd + sd.Nout > 0:
            fo.write(f"\nsd_{suffix[s]} = {sd.sd}, N_{suffix[s]} = {sd.Nsd}, N_out_{suffix[s]} = {sd.Nout}")
        s += 1


def out_file_name_generation(full_out_folder, out_file_name, J, mol_name):
    out_file_name_split = out_file_name.rsplit('.', 1)
    if len(out_file_name_split) == 2:
        out_file_name_out = ".".join([out_file_name_split[0], mol_name, "J" + J, out_file_name_split[1]])
    else:
        out_file_name_out = ".".join([out_file_name_split[0], mol_name, "J" + J])
    out_file_name_full = os.path.join(full_out_folder, out_file_name_out)
    return out_file_name_full


if __name__ == '__main__':
    input_folder = "input"
    output_folder = "output"
    out_file_for_sds = "out_file_all_sd"
    file_calc_name_prefix = "fort.14"

    mode = sys.argv[1]
    mol_name = sys.argv[2]
    file_exp_name = sys.argv[3]
    out_file_exp_name = sys.argv[4]
    out_file_comp_name = sys.argv[5]
    E_zero = float(sys.argv[6])
    E_tr_l = float(sys.argv[7])
    E_tr_h = float(sys.argv[8])
    eps = float(sys.argv[9])
    make_comp_files = sys.argv[10]

    list_inp_files = []
    full_inp_folder = os.path.join(os.getcwd(), input_folder, mol_name)
    if os.path.isdir(full_inp_folder):
        list_inp_files = os.listdir(full_inp_folder)
    else:
        print("Input folder doesn't exist")
        exit(1)

    full_out_folder = os.path.join(os.getcwd(), output_folder, mol_name)
    is_out_dir_exists = os.path.isdir(full_out_folder)
    if not is_out_dir_exists:
        os.makedirs(full_out_folder)

    forts = []
    if mode == "all":
        for file_name in list_inp_files:
            if file_name.startswith(file_calc_name_prefix):
                forts.append(os.path.join(full_inp_folder, file_name))
    else:
        forts.append(os.path.join(full_inp_folder, mode))

    with open(os.path.join(full_out_folder, out_file_for_sds + "_eps=" + str(eps) + "_El=" + str(E_tr_l) + "_Eh=" + str(E_tr_h)+ ".txt"), 'a') as f_sd:
        forts.sort()
        for fort in forts:
            format_calc = formats.CalcFormat(fort, E_zero)
            states_calc = format_calc.parse_file()

            J_list = [str(states_calc[0].J)]

            file_exp_name_full = os.path.join(full_inp_folder, file_exp_name)

            if mol_name.upper() == "H2-16O":
                # for H2-16O
                # 0 ~ v1  1 ~ v2  2 ~ v3  3 ~ J  4 ~ Ka  5 ~ Kc  6 ~ E
                format_exp = formats.ExpFormatH216O(file_exp_name_full, J_list, 3)
            elif mol_name.upper() == "N2O":
                # for N2O
                # 0 ~ q1  1 ~ q2  2 ~ q3  3 ~ sym  4 ~ J  5 ~ E
                format_exp = formats.ExpFormatN2O(file_exp_name_full, J_list, 4)
            else:
                print("Only H2-16O and N2O molecules are supported by now")
                exit(1)

            states_exp = format_exp.parse_file()
            comp_list = do_comparison(states_calc, states_exp, eps)

            out_file_comp_name_full = out_file_name_generation(full_out_folder, out_file_comp_name, J_list[0], mol_name)
            out_file_exp_name_full = out_file_name_generation(full_out_folder, out_file_exp_name, J_list[0], mol_name)

            if make_comp_files == 'True':
                comp_list.write_to_file(out_file_comp_name_full)
                states_exp.write_to_file(out_file_exp_name_full)

            if os.path.exists(out_file_comp_name_full):
                comp_states = comp_list.parse_file(out_file_comp_name_full)
            else:
                print("Can't find file with the comparison! Change the 'make_comp_files' flag to True")
                exit(1)

            out_file_comp_name_split = out_file_comp_name_full.rsplit('.', 1)
            if len(out_file_comp_name_split) == 2:
                out_file_comp_name_full_sd = out_file_comp_name_split[0] + "+sd." + out_file_comp_name_split[1]
            else:
                out_file_comp_name_full_sd = out_file_comp_name_split[0] + "+sd"

            shutil.copyfile(out_file_comp_name_full, out_file_comp_name_full_sd)
            sd_l, sd_m, sd_h = sd_calculation(comp_states, E_tr_l, E_tr_h)
            with open(out_file_comp_name_full_sd, 'a') as f:
                write_sds_to_file(f, [sd_l, sd_m, sd_h])


            f_sd.write(f"\n\nJ = {J_list[0]}")
            write_sds_to_file(f_sd, [sd_l, sd_m, sd_h])


