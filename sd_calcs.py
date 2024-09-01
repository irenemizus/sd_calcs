import sys
from marvel_n2o import get_eners


def parse_fort14(file_name, E_zero):
    au_to_cm = 219474.624
    list_syms = []
    sym = 0
    with open(file_name, 'r') as f:
        line = f.readline()
        while line != '':
            try:
                J, n_levs = int(line.strip().split()[0]), int(line.strip().split()[5])
                sym += 1
                list_ens = []
                r = n_levs // 4 + 1 if n_levs % 4 != 0 else n_levs // 4
                for l in range(r):
                    line = f.readline().replace('D', 'E')
                    list_ens.extend(list(map(lambda x: float(x) * au_to_cm - E_zero, line.strip().split())))
                list_syms.append((sym, list_ens))
            except ValueError:
                print("Something went wrong")

            line = f.readline()

    return J, list_syms


def collect_sym_types(ens_list):
    sym_types = []
    for l in ens_list:
        if l[1] not in sym_types: sym_types.append(l[1])

    return sym_types


if __name__ == '__main__':
    E_zero = float(sys.argv[2])

    file_name = sys.argv[1]
    J, list_syms = parse_fort14(file_name, E_zero)

    exp_ens_list = get_eners(sys.argv[3], [ str(J) ])
    exp_sym_types = collect_sym_types(exp_ens_list)
    list_syms_exp = []

    for s in exp_sym_types:
        list_syms_exp_tmp = []
        for el in exp_ens_list:
            if el[1] == s:
                list_syms_exp_tmp.append(el[2])
        list_syms_exp.append((int(s), list_syms_exp_tmp))



    pass

