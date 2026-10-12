import argparse
import copy
import math
import os.path
import shutil

import formats
import states


def exp_number_det(states_calc: states.States, state_exp, eps, obs_calc_comp):
    # A function for determining a serial number of a calculated energy, which corresponds to the given experimental state
    sign = lambda x: math.copysign(1, x)
    status = states.Status.NOT_FOUND
    list_diff = []
    indx_next = 0
    diff_oc_first = math.nan
    sc_next = None

    if obs_calc_comp:   # an ordinary case of experiment vs calculated data comparison
        E_exp = state_exp.E
    else:               # a case of comparison between two calculated datasets in terms of their obs-calc deviations
        E_exp = state_exp.E_exp

    for sc in states_calc:
        indx_next += 1
        if state_exp.J == sc.J and state_exp.sym == sc.sym:
            sc_next = sc
            diff_oc_first = E_exp - sc_next.E
            list_diff.append((sc_next, diff_oc_first))  # the first candidate
            break

    for i in range(indx_next, len(states_calc)):
        sc_next = states_calc[i]
        if state_exp.J == sc_next.J and state_exp.sym == sc_next.sym:
            diff_oc_next = E_exp - sc_next.E
            list_diff.append((sc_next, diff_oc_next))  # other candidates

            # the one, near which the sign of discrepancy changes, or the one, after which the absolute value of deviation increases, fits best
            if sign(diff_oc_first) != sign(diff_oc_next) or abs(diff_oc_first) <= abs(diff_oc_next):
                break

    if sc_next is None:
        # No calculated state (left) for this (J, sym): this experimental level has no calculated
        # counterpart, which happens when the calculated set is smaller than the experimental one for
        # the given J (a truncated/incomplete calculated dataset). Mark it NOT_FOUND and continue
        # instead of aborting the whole run.
        print(f"WARNING: no calculated state for J = {state_exp.J}, symmetry "
              f"{states.SymType(state_exp.sym).name}, E_exp = {E_exp:.6f} 1/cm (the calculated set is "
              f"smaller than the experimental one for this J); marked as NOT_FOUND and skipped.")
        return None, states.ComparedState(0.0, E_exp, state_exp.J, state_exp.sym, 0, 0.0, 0.0, status, state_exp.qn,
                                          is_hp=state_exp.is_hp, prec_weight=state_exp.prec_weight)

    if sc_next.N == len(states_calc):
        # nothing was found
        return None, states.ComparedState(0.0, E_exp, state_exp.J, state_exp.sym, 0, 0.0, 0.0, status, state_exp.qn,
                                          is_hp=state_exp.is_hp, prec_weight=state_exp.prec_weight)

    if list_diff:
        # the state we're searching for is the one, which discrepancy is the least
        if len(list_diff) == 1:
            # only one candidate remains for this (J, sym); it is the best fit by default
            sc_fit, diff_oc_fit = list_diff[0]
        else:
            (sc_fit, diff_oc_fit) = list_diff[-1] if abs(list_diff[-1][1]) <= abs(list_diff[-2][1]) else list_diff[-2]
        if abs(diff_oc_fit) <= eps:
            # FOUND; in '--rel_weights' mode the weight carries the level reliability grade
            return sc_fit, states.ComparedState(sc_fit.E, E_exp, state_exp.J, state_exp.sym, sc_fit.N, state_exp.prec_weight, diff_oc_fit, states.Status.FOUND, state_exp.qn,
                                                is_hp=state_exp.is_hp, prec_weight=state_exp.prec_weight)
        elif abs(diff_oc_fit) > eps:
            # OUTLIER
            return sc_fit, states.ComparedState(sc_fit.E, E_exp, state_exp.J, state_exp.sym, sc_fit.N, 0.0, diff_oc_fit, states.Status.OUTLIER, state_exp.qn,
                                                is_hp=state_exp.is_hp, prec_weight=state_exp.prec_weight)


def do_comparison(states_calc, states_exp, eps, obs_calc_comp=True):
    states_calc_cpy = copy.deepcopy(states_calc)
    comp_states = []
    # iterating by experimental states
    for se in states_exp:
        # searching for the best fit among the calculated energies
        sc_fit, compared_state = exp_number_det(states_calc_cpy, se, eps, obs_calc_comp)
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
                sd_l.sd += cs.E_diff**2 * cs.w
            else:
                sd_l.Nout += 1
                sd_l.sd_out += cs.E_diff**2
        elif cs.E_exp >= E_tr_h:
            if cs.w > 0:
                sd_h.Nsd += 1
                sd_h.sd += cs.E_diff**2 * cs.w
            else:
                sd_h.Nout += 1
                sd_h.sd_out += cs.E_diff**2
        else:
            if cs.w > 0:
                sd_m.Nsd += 1
                sd_m.sd += cs.E_diff**2 * cs.w
            else:
                sd_m.Nout += 1
                sd_m.sd_out += cs.E_diff**2

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


def assign_precision_weights(states_exp, dE_th):
    # Translates each level uncertainty to a per-transition experimental error dE = unc * sqrt(#lines),
    # splits the levels into high-precision (dE <= dE_th) and low-precision (dE > dE_th) groups, and assigns
    # a reliability weight: 1.0 for HP levels, dE_max/dE for LP levels (dE_max = largest dE among HP levels).
    for st in states_exp:
        st.dE = st.unc * math.sqrt(st.nlines) if st.nlines > 0 else 0.0

    hp_levs = [st for st in states_exp if st.dE <= dE_th]
    dE_max = max((st.dE for st in hp_levs), default=0.0)
    if dE_max <= 0.0:
        dE_max = dE_th

    for st in states_exp:
        if st.dE <= dE_th:
            st.is_hp = True
            st.prec_weight = 1.0
        else:
            st.is_hp = False
            st.prec_weight = dE_max / st.dE

    n_hp = len(hp_levs)
    n_lp = len(states_exp) - n_hp
    print(f"--rel_weights: dE_th = {dE_th} 1/cm, dE_max = {dE_max:.6g} 1/cm, "
          f"high-precision = {n_hp}, low-precision = {n_lp}")
    return dE_max


def sd_calculation_rel(comp_states, E_tr_l, E_tr_h):
    # Two-group sd: for each energy region (l/m/h) and each precision group (hp/lp), a StandDev is computed.
    # FOUND states are treated as weight 1.00 (unweighted); OUTLIER states go to sd_out as in sd_calculation.
    # Returns [hp_l, lp_l, hp_m, lp_m, hp_h, lp_h].
    sd = {(g, r): StandDev() for g in ('hp', 'lp') for r in ('l', 'm', 'h')}

    for cs in comp_states:
        if cs.E_exp <= E_tr_l:
            region = 'l'
        elif cs.E_exp >= E_tr_h:
            region = 'h'
        else:
            region = 'm'
        group = 'hp' if cs.is_hp else 'lp'
        obj = sd[(group, region)]
        if cs.w > 0:
            obj.Nsd += 1
            obj.sd += cs.E_diff**2          # all non-zero weights are taken as 1.00
        else:
            obj.Nout += 1
            obj.sd_out += cs.E_diff**2

    out = []
    for region in ('l', 'm', 'h'):
        for group in ('hp', 'lp'):
            obj = sd[(group, region)]
            if obj.Nsd > 0:
                obj.sd = math.sqrt(obj.sd / obj.Nsd)
            if obj.Nout > 0:
                obj.sd_out = math.sqrt(obj.sd_out / obj.Nout)
            out.append(obj)

    return out


def calc_sd(comp_states, E_tr_l, E_tr_h, use_rel):
    if use_rel:
        return sd_calculation_rel(comp_states, E_tr_l, E_tr_h)
    return list(sd_calculation(comp_states, E_tr_l, E_tr_h))


def _write_one_sd(fo, label, sd):
    if sd.Nsd + sd.Nout > 0:
        Nout_rel = sd.Nout / (sd.Nsd + sd.Nout) * 100
        fo.write(f"\nsd_{label} = {sd.sd:.4f}, N_{label} = {sd.Nsd}, sd_{label}_out = {sd.sd_out:.4f}, "
                 f"N_out_{label} = {sd.Nout}, N_out_rel_{label} = {Nout_rel:.2f}%")


def write_sds_to_file(fo, sd_list):
    suffix = ["l", "m", "h"]
    s = 0
    for sd in sd_list:
        _write_one_sd(fo, suffix[s], sd)
        s += 1


def write_sds_to_file_rel(fo, sd_list):
    labels = ["hp_l", "lp_l", "hp_m", "lp_m", "hp_h", "lp_h"]
    for label, sd in zip(labels, sd_list):
        _write_one_sd(fo, label, sd)


def write_sds(fo, sd_list, use_rel):
    if use_rel:
        write_sds_to_file_rel(fo, sd_list)
    else:
        write_sds_to_file(fo, sd_list)


def outliers_by_N_out_rel(comp_states, N_out_rel):
    comp_states_out = []

    def w_zero(x):
        x.w, x.status = 0.0, states.Status.OUTLIER
        return x

    def w_one(x):
        # keep the reliability grade (1.0 unless '--rel_weights' set a lower value)
        x.w, x.status = x.prec_weight, states.Status.FOUND
        return x

    comp_states_sort_by_diff = sorted(comp_states, key=lambda x: abs(x.E_diff))
    N_out = math.floor(len(comp_states) * N_out_rel / 100.0)
    N_found = len(comp_states) - N_out
    found_body = list(map(w_one, comp_states_sort_by_diff[:N_found]))   # All the levels with E_diff <= eps_max become FOUND
    comp_states_out.extend(found_body)
    eps_max = comp_states_out[-1].E_diff
    outliers_tail = list(map(w_zero, comp_states_sort_by_diff[N_found:])) # The N_out_rel % remaining levels become OUTLIERs
    comp_states_out.extend(outliers_tail)
    return sorted(comp_states_out, key=lambda x: (x.J, x.sym, x.E_exp)), eps_max


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
        - files with calculated data: fort.14-like, ExoMol states-like
        - files with observed data: MARVEL-like 
        - output formats: Yurchenko's fit input file
        - output labeling formats: A/B symmetries (A1, A2, B1, B2), 6 quantum numbers for H2-16O, 4 quantum numbers for N2O  
    """

    # Usually non-changeable parameters
    au_to_cm = 219474.624                   # A coefficient for transferring an energy value from E_h to 1 / cm
    input_folder = "input"                  # A folder containing input data (fort.14 files with calculated energies and a file with experimental data)
    output_folder = "output"                # A folder containing output data
    out_file_for_sds = "out_file_all_sd"    # A name for an output file containing only sd values for all the J values and energy ranges of interest (without extension)
    file_calc_name_prefix = "fort.14"       # A starting part of the names for the fort.14-like input files with calculated energies

    # Command line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", type=str,
                        help="specifies if the script should process all the 'fort.14' files ('all' mode) "
                             "or only the specified one (including a single ExoMol-format states file)",
                        default='all')
    parser.add_argument("--mol_name", type=str, help="name of the molecule", choices=['H2-16O', 'N2O', 'N2O-556'],
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
    parser.add_argument("--Jmax", type=int,
                        help="maximum value of J for comparison with ExoMol states-like calculated data", default=100)

    parser.add_argument("--make_comp_files", type=str, choices=['True', 'False'],
                        help="defines if the output 'out_file_exp_name' and 'out_file_comp_name' files should be generated ('True' mode), "
                             "or the existed ones should be taken from the 'output_folder/mol_name' directory ('False' mode)", default='True')

    parser.add_argument("--N_out_rel", type=float,
                        help="a relative number of worst-predicted states, which should be considered as outliers for a final all-J sd calculation, in percentage", default=0.0)

    parser.add_argument("--rel_weights", action='store_true',
                        help="assign relative weights to the levels based on the experimental uncertainties "
                             "(requires the input file with experimental data to carry the 'unc' and '#lines' columns, "
                             "as in the ExoMol/MARVEL and RITZ formats). High-precision levels get weight 1.0, "
                             "low-precision levels get weight dE_max/dE. The sd is then reported separately for the "
                             "high- and low-precision groups. Requires --make_comp_files True.")
    parser.add_argument("--dE_th", type=float,
                        help="per-transition experimental error threshold (1/cm) separating high-precision (dE <= dE_th) "
                             "from low-precision (dE > dE_th) levels; only used with --rel_weights", default=0.001)

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
    Jmax = args.Jmax
    N_out_rel = args.N_out_rel
    make_comp_files = args.make_comp_files
    rel_weights = args.rel_weights
    dE_th = args.dE_th
    # finer weight column resolution is needed to keep the fractional reliability grades
    w_spec = '7.4f' if rel_weights else '5.2f'

    if rel_weights and make_comp_files != 'True':
        print("ERROR: '--rel_weights' requires '--make_comp_files True' (the high/low-precision split must be "
              "known at the sd-calculation step).")
        exit(1)

    # Collecting input files
    list_inp_files = []
    full_inp_folder = os.path.join(os.getcwd(), input_folder, mol_name)
    if os.path.isdir(full_inp_folder):
        list_inp_files = os.listdir(full_inp_folder)
    else:
        print(f"Input folder '{full_inp_folder}' doesn't exist")
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
        comp_states_allJ = []
        rel_warned = False
        for fort in forts:
            if not mode.startswith("fort.14") and len(forts) == 1 and not fort.split("/")[-1].startswith("fort.14"):
                J_list = [str(J) for J in range(Jmax + 1)]
                # For ExoMol (e.g., POKAZATEL) data ----------------
                # Obtaining calculated states
                exomol_format = formats.ExoMolFormatH216O(fort, J_list, 3)
                states_calc = exomol_format.parse_file()
                # ------------------------------------
            else:
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
            elif mol_name.upper() == "N2O-556":
                # for N2O-556 (RITZ format); J is the first column
                format_exp = formats.RITZFormatN2O556(file_exp_name_full, J_list, 0)
            else:
                print("Only H2-16O, N2O, and N2O-556 molecules are supported by now")
                exit(1)

            states_exp = format_exp.parse_file()

            # '--rel_weights': only active if the experimental data actually carries uncertainties
            has_unc = any(st.nlines > 0 or st.unc > 0 for st in states_exp)
            use_rel = rel_weights and has_unc
            if rel_weights and not has_unc and not rel_warned:
                print("WARNING: '--rel_weights' passed but the experimental file has no 'unc'/'#lines' columns; "
                      "falling back to standard weighting.")
                rel_warned = True
            if use_rel:
                assign_precision_weights(states_exp, dE_th)

            # Generating full output file names
            if not mode.startswith("fort.14") and len(forts) == 1:
                J_suffix = f"all_eps={eps}"
            else:
                J_suffix = J_list[0] + f"_eps={eps}"

            out_file_comp_name_full = out_file_name_generation(full_out_folder, out_file_comp_name, J_suffix, mol_name)
            out_file_exp_name_full = out_file_name_generation(full_out_folder, out_file_exp_name, J_suffix, mol_name)
            out_file_comp_name_full_allJ = out_file_name_generation(full_out_folder, out_file_comp_name, f'all_eps={eps}', mol_name)

            comp_list = states.ComparisonList([])
            # Generating output files with comparison results (for 'True' mode)
            if make_comp_files == 'True':
                # Obtaining a list of compared states
                comp_list_init = do_comparison(states_calc, states_exp, eps)
                comp_list = states.ComparisonList([])
                comp_list.add_item(comp_list_init[0])
                for i in range(len(comp_list_init) - 1):
                    prev_line = comp_list_init[i]
                    cur_line = comp_list_init[i + 1]
                    # N is assigned per (J, symmetry) block and resets when the block changes, so the
                    # monotonic-increasing check only applies within the same block
                    if (prev_line.J != cur_line.J) or (prev_line.sym != cur_line.sym) or (cur_line.N > prev_line.N):
                        comp_list.add_item(cur_line)
                    else:
                        print(f"WARNING: There is a line with J = {cur_line.J}, symmetry {cur_line.sym},"
                              f" experimental energy {cur_line.E_exp} 1/cm "
                              f"and the corresponding calculated energy {cur_line.E_calc} 1/cm in the "
                              f"comparison list, whose number N = {cur_line.N} is less than "
                              f"the one of the previous line N_prev = {prev_line.N}. "
                              f"Skipping it. Please check your results. "
                              f"Maybe, the amount of the calculated data is insufficient.")

                comp_list.write_to_file(out_file_comp_name_full, w_spec=w_spec)
                if len(forts) > 1:
                    comp_list.write_to_file(out_file_comp_name_full_allJ, mode='a', w_spec=w_spec)
                states_exp.write_to_file(out_file_exp_name_full, w_spec=w_spec)

            """
                Here there is a possibility to change somehow the content of the pre-generated file with comparison results.
                For example, to change weights for some of the compared states before calculating sd for them.
                That's why 'make_comp_files' option is needed,  
                Besides, in the 'False' mode all the states with zero weights (w = 0.0) will be written to a separate file 
                with '+sd_valid' suffix in its name. 
            """

            # Parsing the pre-generated output file with comparison results
            if use_rel:
                # use the in-memory comparison list so the high/low-precision (is_hp) grouping is preserved
                comp_states = list(comp_list)
                comp_states_allJ.extend(comp_states)
            elif os.path.exists(out_file_comp_name_full):
                comp_states = comp_list.parse_file(out_file_comp_name_full, format_exp.name)
                comp_states_allJ.extend(comp_states)
            else:
                print("Can't find file with the comparison! Change the 'make_comp_files' flag to True")
                exit(1)

            # Generating a name with '+sd' suffix for a new output file, which will contain both the comparison and sd values
            out_file_comp_name_split = out_file_comp_name_full.rsplit('.', 1)
            if len(out_file_comp_name_split) == 2:  # if out_file_comp_name_full contains extension
                out_file_comp_name_full_sd = out_file_comp_name_split[0] + f"+sd." + out_file_comp_name_split[1]
            else:
                out_file_comp_name_full_sd = out_file_comp_name_split[0] + f"+sd"

            # Copying comparison results to the new file
            shutil.copyfile(out_file_comp_name_full, out_file_comp_name_full_sd)

            # Calculating sds and writing them at the end of the new file
            sd_list = calc_sd(comp_states, E_tr_l, E_tr_h, use_rel)
            with open(out_file_comp_name_full_sd, 'a') as f:
                write_sds(f, sd_list, use_rel)

            # Generating a name with '+sd_valid' suffix for another output file, which will contain comparison and sds only for the zero-weighed states
            if make_comp_files == 'False':
                print("NOTE: changing the '--eps' parameter requires the '--make_comp_files' parameter to be True! "
                      "Recalculation from the very beginning is required for every new eps value!")
                if len(out_file_comp_name_split) == 2:  # if out_file_comp_name_full contains extension
                    out_file_comp_name_full_sd_val = out_file_comp_name_split[0] + f"+sd_eps={eps}_valid." + \
                                                     out_file_comp_name_split[1]
                else:
                    out_file_comp_name_full_sd_val = out_file_comp_name_split[0] + f"+sd_eps={eps}_valid"

                # Copying comparison results for the zero-weighed states to the new file
                comp_states_val = comp_list.write_to_file(out_file_comp_name_full_sd_val, filter_by=0.0)

                # Calculating sds for the zero-weighed states and writing them at the end of the new file
                # (this path is non-rel only: '--rel_weights' requires '--make_comp_files True')
                sd_val_list = calc_sd(comp_states_val, E_tr_l, E_tr_h, use_rel=False)
                with open(out_file_comp_name_full_sd_val, 'a') as f_val:
                    write_sds(f_val, sd_val_list, use_rel=False)

            # Collecting sd values for all the J values of interest in one file
            if len(forts) > 1:
                f_sd.write(f"\n\nJ = {J_list[0]}")
                write_sds(f_sd, sd_list, use_rel)

        # sd calculation for all J values
        sd_list_allJ = calc_sd(comp_states_allJ, E_tr_l, E_tr_h, use_rel)
        f_sd.write(f"\n\nall Js:")
        write_sds(f_sd, sd_list_allJ, use_rel)

        if N_out_rel > 0.0:     # if the worst N_out_rel % levels need to be considered OUTLIERs
            # eps_max corresponds to N_out_rel % OUTLIERs
            comp_states_allJ_wo_N_rel_lvls, eps_max = outliers_by_N_out_rel(comp_states_allJ, N_out_rel)

            # writing the new set of levels to a file
            out_file_comp_name_full_allJ_wo_Nrel_lvls = out_file_name_generation(full_out_folder, out_file_comp_name,
                                                                                 f'all+sd_eps={eps}_wo_Nrel={N_out_rel}%_lvls', mol_name)
            comp_list_allJ_wo_N_rel_lvls = states.ComparisonList(comp_states_allJ_wo_N_rel_lvls)
            comp_list_allJ_wo_N_rel_lvls.write_to_file(out_file_comp_name_full_allJ_wo_Nrel_lvls, w_spec=w_spec)

            # calculating new sd for all J values
            sd_list_wo_Nrel_lvls = calc_sd(comp_states_allJ_wo_N_rel_lvls, E_tr_l, E_tr_h, use_rel)
            # appending the new sd values to the file
            with open(out_file_comp_name_full_allJ_wo_Nrel_lvls, 'a') as f_Nout:
                f_Nout.write(f"\n\nall Js, eps_max = {eps_max:5.3f}:")
                write_sds(f_Nout, sd_list_wo_Nrel_lvls, use_rel)

            # writing the new set of outliers to a file
            out_file_comp_name_full_allJ_wo_Nrel_lvls_val = out_file_name_generation(full_out_folder, out_file_comp_name,
                                                                                  f'all+sd_eps={eps}_valid_wo_Nrel={N_out_rel}%_lvls', mol_name)
            comp_out_states_val = comp_list_allJ_wo_N_rel_lvls.write_to_file(out_file_comp_name_full_allJ_wo_Nrel_lvls_val, filter_by=0.0, w_spec=w_spec)

            # calculating new sd for all J values
            sd_list_wo_Nrel_lvls_val = calc_sd(comp_out_states_val, E_tr_l, E_tr_h, use_rel)
            # appending the new sd values to the file
            with open(out_file_comp_name_full_allJ_wo_Nrel_lvls_val, 'a') as f_Nout_val:
                f_Nout_val.write(f"\n\nall Js, eps_max = {eps_max:5.3f}:")
                write_sds(f_Nout_val, sd_list_wo_Nrel_lvls_val, use_rel)




