import math

from . import states


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


def assign_precision_weights(states_exp, dE_th, hp_weight=1.0):
    # Translates each level uncertainty to a per-transition experimental error dE = unc * sqrt(#lines),
    # splits the levels into high-precision (dE <= dE_th) and low-precision (dE > dE_th) groups, and assigns
    # a reliability weight: hp_weight for HP levels, dE_max/dE for LP levels (dE_max = largest dE among HP levels).
    for st in states_exp:
        st.dE = st.unc * math.sqrt(st.nlines) if st.nlines > 0 else 0.0

    hp_levs = [st for st in states_exp if st.dE <= dE_th]
    dE_max = max((st.dE for st in hp_levs), default=0.0)
    if dE_max <= 0.0:
        dE_max = dE_th

    for st in states_exp:
        if st.dE <= dE_th:
            st.is_hp = True
            st.prec_weight = hp_weight
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
