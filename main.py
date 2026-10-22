"""
A script for finding related pairs of calculated and observed states for different molecules and obtaining
deviations and sd values for them.

By now, the following is supported:
- molecules: H2-16O, N2O, N2O-556
- files with calculated data: fort.14-like, ExoMol states-like, HITRAN
- files with observed data: MARVEL-like, RITZ-like (for `N2O-556`)
- output formats: Yurchenko's fit input file, comparison files, sd summary files
- output labeling formats:
  - symmetries `A1, A2, B1, B2` (plus `UKN` for undetermined)
  - 6 quantum numbers for `H2-16O` (`v1, v2, v3, J, Ka, Kc`)
  - 4 quantum numbers for `N2O` (`v1, v2, J, l`)
  - 5 quantum numbers for `N2O-556` (`v1, v2, v3, J, l`)
"""

import argparse
import os.path
import shutil

from sdcalcs import formats
from sdcalcs import states
from sdcalcs.comparison import do_comparison
from sdcalcs.sd import (assign_precision_weights, calc_sd, write_sds, outliers_by_N_out_rel)
from sdcalcs.utils import out_file_name_generation


if __name__ == '__main__':

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
    parser.add_argument("--hp_weight", type=float,
                        help="weight assigned to high-precision levels; only used with --rel_weights", default=1.0)
    parser.add_argument("--E_cut", type=float,
                        help="cut the ens fit-input files at this energy (1/cm): levels with E > E_cut are omitted",
                        default=None)
    parser.add_argument("--E_fit_max", type=float,
                        help="maximum energy (1/cm) for the fit: levels with E > E_fit_max in the "
                             "ens_Yur_format.* files get zero weight; all other output files are unaffected",
                        default=None)

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
    hp_weight = args.hp_weight
    E_cut = args.E_cut
    E_fit_max = args.E_fit_max
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
        suboptimal_all = []
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
                assign_precision_weights(states_exp, dE_th, hp_weight)

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
                        suboptimal_all.append(cur_line)
                        print(f"WARNING: suboptimal match: J = {cur_line.J}, symmetry {cur_line.sym},"
                              f" E_exp = {cur_line.E_exp:.6f} 1/cm, E_calc = {cur_line.E_calc:.6f} 1/cm, "
                              f"N = {cur_line.N} < N_prev = {prev_line.N}. "
                              f"Written to suboptimal_matches.txt.")

                comp_list.write_to_file(out_file_comp_name_full, w_spec=w_spec)
                if len(forts) > 1:
                    comp_list.write_to_file(out_file_comp_name_full_allJ, mode='a', w_spec=w_spec)
                # The Yurchenko fit-input (ens) file must contain exactly the same set of levels as the
                # comparison file written just above: only the matched states that passed the monotonic-
                # increasing N check. states_exp still holds every experimental state (including the
                # NOT_FOUND ones with N = 0 and the non-monotonic ones), so it is truncated down to the
                # levels actually present in comp_list. Each level's quantum-number object is unique and
                # is shared between its State and ComparedState, so membership is tested by that object.
                kept_qn = {cs.qn for cs in comp_list}
                ens_states = [se for se in states_exp if se.qn in kept_qn]
                if E_cut is not None:
                    ens_states = [se for se in ens_states if se.E <= E_cut]
                if E_fit_max is not None:
                    for se in ens_states:
                        if se.E > E_fit_max:
                            se.w = 0.0
                states.States(ens_states).write_to_file(out_file_exp_name_full, w_spec=w_spec)

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

        # Writing suboptimal (non-monotonic N) matches to a separate file
        if suboptimal_all:
            suboptimal_file = os.path.join(full_out_folder, 'suboptimal_matches.txt')
            states.ComparisonList(suboptimal_all).write_to_file(suboptimal_file, w_spec=w_spec)
            print(f"\n{len(suboptimal_all)} suboptimal match(es) written to {suboptimal_file}")

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
