import argparse
import copy
import math
import os.path
import shutil

import formats
import states


def exp_number_det(states_calc: states.States, state_exp, eps):
    # A function for determining a serial number of a calculated energy, which corresponds to the given experimental state
    sign = lambda x: math.copysign(1, x)
    status = states.Status.NOT_FOUND
    list_diff = []
    indx_next = 0
    diff_oc_first = math.nan
    sc_next = None

    for sc in states_calc:
        indx_next += 1
        if state_exp.J == sc.J and state_exp.sym == sc.sym:
            sc_next = sc
            diff_oc_first = state_exp.E - sc_next.E
            list_diff.append((sc_next, diff_oc_first))  # the first candidate
            break

    for i in range(indx_next, len(states_calc)):
        sc_next = states_calc[i]
        if state_exp.J == sc_next.J and state_exp.sym == sc_next.sym:
            diff_oc_next = state_exp.E - sc_next.E
            list_diff.append((sc_next, diff_oc_next))  # other candidates

            # the one, near which the sign of discrepancy changes, or the one, after which the absolute value of deviation increases, fits best
            if sign(diff_oc_first) != sign(diff_oc_next) or abs(diff_oc_first) <= abs(diff_oc_next):
                break

    if sc_next is None:
        print("Something wrong is happening! Please check the file with experimental data!")
        exit(1)

    if sc_next.N == len(states_calc):
        # nothing was found
        return None, states.ComparedState(0.0, state_exp.E, state_exp.J, state_exp.sym, 0, 0.0, 0.0, status, state_exp.qn)

    if list_diff:
        # the state we're searching for is the one, which discrepancy is the least
        (sc_fit, diff_oc_fit) = list_diff[-1] if abs(list_diff[-1][1]) <= abs(list_diff[-2][1]) else list_diff[-2]
        if abs(diff_oc_fit) <= eps:
            # FOUND
            return sc_fit, states.ComparedState(sc_fit.E, state_exp.E, state_exp.J, state_exp.sym, sc_fit.N, 1.0, diff_oc_fit, states.Status.FOUND, state_exp.qn)
        elif abs(diff_oc_fit) > eps:
            # OUTLIER
            return sc_fit, states.ComparedState(sc_fit.E, state_exp.E, state_exp.J, state_exp.sym, sc_fit.N, 0.0, diff_oc_fit, states.Status.OUTLIER, state_exp.qn)


def do_comparison(states_calc, states_exp, eps):
    states_calc_cpy = copy.deepcopy(states_calc)
    comp_states = []
    # iterating by experimental states
    for se in states_exp:
        # searching for the best fit among the calculated energies
        sc_fit, compared_state = exp_number_det(states_calc_cpy, se, eps)
        if sc_fit:
            states_calc_cpy.remove_item(sc_fit)  # removing the found calculated energy
            se.N = compared_state.N              # saving the serial number of the found energy...
            se.w = compared_state.w              # ...and its status (FOUND/OUTLIER) to the current experimental state
            comp_states.append(compared_state)
    return states.ComparisonList(comp_states)


class StandDev:
    def __init__(self, Nsd=0, Nout=0, sd=0.0, sd_out=0.0):
        self.sd = sd            # standard deviation for the FOUND states
        self.sd_out = sd_out    # standard deviation for the OUTLIER states
        self.Nsd = Nsd          # number of FOUND states
        self.Nout = Nout        # number of OUTLIER states


def sd_calculation(comp_states, E_tr_l, E_tr_h):
    sd_l = StandDev()           # StandDev object for low-energy states (E <= E_tr_l)
    sd_m = StandDev()           # StandDev object for middle-energy states (E_tr_l < E < E_tr_h)
    sd_h = StandDev()           # StandDev object for high-energy states (E >= E_tr_h)

    for cs in comp_states:
        if cs.E_exp <= E_tr_l:
            if cs.w > 0:
                sd_l.Nsd += 1
                sd_l.sd += cs.E_diff ** 2 * cs.w
            else:
                sd_l.Nout += 1
                sd_l.sd_out += cs.E_diff ** 2
        elif cs.E_exp >= E_tr_h:
            if cs.w > 0:
                sd_h.Nsd += 1
                sd_h.sd += cs.E_diff**2 * cs.w
            else:
                sd_h.Nout += 1
                sd_h.sd_out += cs.E_diff ** 2
        else:
            if cs.w > 0:
                sd_m.Nsd += 1
                sd_m.sd += cs.E_diff**2 * cs.w
            else:
                sd_m.Nout += 1
                sd_m.sd_out += cs.E_diff ** 2

    if sd_h.Nsd > 0:
        sd_h.sd = math.sqrt(sd_h.sd / sd_h.Nsd)
    if sd_m.Nsd > 0:
        sd_m.sd = math.sqrt(sd_m.sd / sd_m.Nsd)
    if sd_l.Nsd > 0:
        sd_l.sd = math.sqrt(sd_l.sd / sd_l.Nsd)

    if sd_h.Nout > 0:
        sd_h.sd_out = math.sqrt(sd_h.sd_out / sd_h.Nout)
    if sd_m.Nout > 0:
        sd_m.sd_out = math.sqrt(sd_m.sd_out / sd_m.Nout)
    if sd_l.Nout > 0:
        sd_l.sd_out = math.sqrt(sd_l.sd_out / sd_l.Nout)

    return sd_l, sd_m, sd_h


def write_sds_to_file(fo, sd_list):
    suffix = ["l", "m", "h"]
    s = 0
    for sd in sd_list:
        if sd.Nsd + sd.Nout > 0:
            Nout_rel = sd.Nout / (sd.Nsd + sd.Nout) * 100
            fo.write(f"\nsd_{suffix[s]} = {sd.sd:.4f}, N_{suffix[s]} = {sd.Nsd}, sd_{suffix[s]}_out = {sd.sd_out:.4f}, N_out_{suffix[s]} = {sd.Nout}, N_out_rel_{suffix[s]} = {Nout_rel:.2f}%")
        s += 1


def out_file_name_generation(full_out_folder, out_file_name, J, mol_name):
    out_file_name_split = out_file_name.rsplit('.', 1)
    if len(out_file_name_split) == 2:  # if out_file_name contains extension
        out_file_name_out = ".".join([out_file_name_split[0], mol_name, "J" + J, out_file_name_split[1]])
    else:
        out_file_name_out = ".".join([out_file_name_split[0], mol_name, "J" + J])
    out_file_name_full = os.path.join(full_out_folder, out_file_name_out)
    return out_file_name_full


if __name__ == '__main__':
    """
        A script for finding related pairs of calculated and observed states for different molecules and obtaining 
        deviations and sd values for them.
        
        By now, the following is supported:
        - molecules: H2-16O, N2O
        - files with calculated data: fort.14-like
        - files with observed data: MARVEL-like 
        - output formats: Yurchenko's fit input file
        - output labeling formats: A/B symmetries (A1, A2, B1, B2), 6 quantum numbers for H2-16O, 4 quantum numbers for N2O  
    """

    # Usually non-changeable parameters
    au_to_cm = 219474.624                   # A coefficient for transferring an energy value from E_h to 1 / cm
    input_folder = "input"                  # A folder containing input data (fort.14 files with calculated energies and a file with experimental data)
    output_folder = "output"                # A folder containing output data
    out_file_for_sds = "out_file_all_sd"    # A name for an output file containing only sd values for all the J values and energy ranges of interest (without extension)
    file_calc_name_prefix = "fort.14"       # A starting part of the names for the input files with calculated energies

    # Command line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", type=str,
                        help="specifies if the script should process all the 'fort.14' files ('all' mode) or only the specified one",
                        default='all')
    parser.add_argument("--mol_name", type=str, help="name of the molecule", choices=['H2-16O', 'N2O'],
                        default='H2-16O')
    parser.add_argument("--file_exp_name", type=str, help="name of the input file with experimental data")
    parser.add_argument("--out_file_exp_name", type=str,
                        help="basic part of the output file name; will be extended with the current parameter values",
                        default="ens_Yur_format.txt")
    parser.add_argument("--out_file_comp_name", type=str,
                        help="basic part of the output, containing the entire comparison results, file name; will be extended with the current parameter values",
                        default="comp.txt")

    parser.add_argument("--E_zero", type=float, help="zero energy value for calculated data")
    parser.add_argument("--E_tr_l", type=float, help="energy threshold value limiting the upper edge for low-energy states", default=15000.0)
    parser.add_argument("--E_tr_h", type=float,
                        help="energy threshold value limiting the upper edge for middle-energy states", default=25000.0)
    parser.add_argument("--eps", type=float,
                        help="calculated energies with |obs-calc| <= eps will be considered as the 'FOUND' ones; the others will be marked as 'OUTLIER's", default=0.2)

    parser.add_argument("--make_comp_files", type=str, choices=['True', 'False'],
                        help="defines if the output 'out_file_exp_name' and 'out_file_comp_name' files should be generated ('True' mode), "
                             "or the existed ones should be taken from the 'output_folder/mol_name' directory ('False' mode)", default='True')

    args = parser.parse_args()

    mode = args.mode
    mol_name = args.mol_name
    file_exp_name = args.file_exp_name
    out_file_exp_name = args.out_file_exp_name
    out_file_comp_name = args.out_file_comp_name

    E_zero = args.E_zero
    E_tr_l = args.E_tr_l
    E_tr_h = args.E_tr_h
    eps = args.eps
    make_comp_files = args.make_comp_files

    # Collecting input files
    list_inp_files = []
    full_inp_folder = os.path.join(os.getcwd(), input_folder, mol_name)
    if os.path.isdir(full_inp_folder):
        list_inp_files = os.listdir(full_inp_folder)
    else:
        print("Input folder doesn't exist")
        exit(1)

    # Creating/checking output directories
    full_out_folder = os.path.join(os.getcwd(), output_folder, mol_name)
    is_out_dir_exists = os.path.isdir(full_out_folder)
    if not is_out_dir_exists:
        os.makedirs(full_out_folder)

    # Making list with input files with calculated energies (for mode='all')...
    forts = []
    if mode == "all":
        for file_name in list_inp_files:
            if file_name.startswith(file_calc_name_prefix):
                forts.append(os.path.join(full_inp_folder, file_name))
    # ...or using the specified by '--mode' file, only
    else:
        forts.append(os.path.join(full_inp_folder, mode))

    # Creating the file for sd values
    with open(os.path.join(full_out_folder, out_file_for_sds + "_eps=" + str(eps) + "_El=" + str(E_tr_l) + "_Eh=" + str(E_tr_h)+ ".txt"), 'w') as f_sd:
        forts.sort()
        for fort in forts:
            # Obtaining calculated states for the given J
            format_calc = formats.CalcFormat(fort, E_zero, au_to_cm)
            states_calc = format_calc.parse_file()

            # Getting the current J value
            J_list = [str(states_calc[0].J)]

            # Getting the full path for the input file with experimental data
            file_exp_name_full = os.path.join(full_inp_folder, file_exp_name)

            # Obtaining experimental states for the given J
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

            # Obtaining a list of compared states
            comp_list = do_comparison(states_calc, states_exp, eps)

            # Generating full output file names
            out_file_comp_name_full = out_file_name_generation(full_out_folder, out_file_comp_name, J_list[0], mol_name)
            out_file_exp_name_full = out_file_name_generation(full_out_folder, out_file_exp_name, J_list[0], mol_name)

            # Generating output files with comparison results (for 'True' mode)
            if make_comp_files == 'True':
                comp_list.write_to_file(out_file_comp_name_full)
                states_exp.write_to_file(out_file_exp_name_full)

            """
                Here there is a possibility to change somehow the content of the pre-generated file with comparison results.
                For example, to change weights for some of the compared states before calculating sd for them.
                That's why 'make_comp_files' option is needed,  
                Besides, in the 'False' mode all the states with zero weights (w = 0.0) will be written to a separate file 
                with '+sd_valid' suffix in its name. 
            """

            # Parsing the pre-generated output file with comparison results
            if os.path.exists(out_file_comp_name_full):
                comp_states = comp_list.parse_file(out_file_comp_name_full)
            else:
                print("Can't find file with the comparison! Change the 'make_comp_files' flag to True")
                exit(1)

            # Generating a name with '+sd' suffix for a new output file, which will contain both the comparison and sd values
            out_file_comp_name_split = out_file_comp_name_full.rsplit('.', 1)
            if len(out_file_comp_name_split) == 2:  # if out_file_comp_name_full contains extension
                out_file_comp_name_full_sd = out_file_comp_name_split[0] + "+sd." + out_file_comp_name_split[1]
            else:
                out_file_comp_name_full_sd = out_file_comp_name_split[0] + "+sd"

            # Copying comparison results to the new file
            shutil.copyfile(out_file_comp_name_full, out_file_comp_name_full_sd)

            # Calculating sds and writing them at the end of the new file
            sd_l, sd_m, sd_h = sd_calculation(comp_states, E_tr_l, E_tr_h)
            with open(out_file_comp_name_full_sd, 'a') as f:
                write_sds_to_file(f, [sd_l, sd_m, sd_h])

            # Generating a name with '+sd_valid' suffix for another output file, which will contain comparison and sds only for the zero-weighed states
            if make_comp_files == 'False':
                if len(out_file_comp_name_split) == 2:  # if out_file_comp_name_full contains extension
                    out_file_comp_name_full_sd_val = out_file_comp_name_split[0] + "+sd_valid." + \
                                                     out_file_comp_name_split[1]
                else:
                    out_file_comp_name_full_sd_val = out_file_comp_name_split[0] + "+sd_valid"

                # Copying comparison results for the zero-weighed states to the new file
                comp_states_val = comp_list.write_to_file(out_file_comp_name_full_sd_val, filter_by=0.0)

                # Calculating sds for the zero-weighed states and writing them at the end of the new file
                sd_val_l, sd_val_m, sd_val_h = sd_calculation(comp_states_val, E_tr_l, E_tr_h)
                with open(out_file_comp_name_full_sd_val, 'a') as f_val:
                    write_sds_to_file(f_val, [sd_val_l, sd_val_m, sd_val_h])

            # Collecting sd values for all the J values of interest in one file
            f_sd.write(f"\n\nJ = {J_list[0]}")
            write_sds_to_file(f_sd, [sd_l, sd_m, sd_h])


