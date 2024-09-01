import copy
import math
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


def sd_calculation(comp_states):
    Nsd = 0
    sd = 0.0
    for cs in comp_states:
        if cs.w > 0:
            Nsd += 1
            sd += cs.E_diff**2 * cs.w
    if Nsd > 0:
        sd = math.sqrt(sd / Nsd)

    return sd, Nsd


if __name__ == '__main__':
    #J_list=['2', '15', '20']
    E_zero = float(sys.argv[3])
    file_calc_name = sys.argv[1]
    file_exp_name = sys.argv[2]
    out_file_exp_name = sys.argv[4]
    out_file_comp_name = sys.argv[5]
    eps = float(sys.argv[6])

    format_calc = formats.CalcFormat(file_calc_name, E_zero)
    states_calc = format_calc.parse_file()

    J_list = [str(states_calc[0].J)]

    # 0 ~ v1  1 ~ v2  2 ~ v3  3 ~ J  4 ~ Ka  5 ~ Kc  6 ~ E
    #format_exp_h216o = formats.ExpFormatH216O(file_exp_name, J_list, 3)
    #states_exp_h216o = format_exp_h216o.parse_file()

    # 0 ~ q1  1 ~ q2  2 ~ q3  3 ~ sym  4 ~ J  5 ~ E
    format_exp_n2o = formats.ExpFormatN2O(file_exp_name, J_list, 4)
    states_exp_n2o = format_exp_n2o.parse_file()

    comp_list = do_comparison(states_calc, states_exp_n2o, eps)
    comp_list.write_to_file(out_file_comp_name)
    states_exp_n2o.write_to_file(out_file_exp_name)

    comp_states = comp_list.parse_file(out_file_comp_name)
    sd, Nsd = sd_calculation(comp_states)
    with open(out_file_comp_name, 'a') as f:
        f.write(f"\nsd = {sd}, N = {Nsd}\n")



