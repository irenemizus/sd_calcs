import copy
import math

from . import states


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
