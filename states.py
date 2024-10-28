from enum import Enum


class State:
    def __init__(self, E, J, sym, N=0, w=1.00, qn=None):
        self.E = E
        self.J = J
        self.sym = sym
        self.N = N
        self.w = w
        self.qn = qn


class QuantNumbers:
    def __init__(self, qn_list):
        self.__qn_list = qn_list

    def __len__(self):
        return len(self.__qn_list)

    def __iter__(self):
        return iter(self.__qn_list)

    def __getitem__(self, item):
        return self.__qn_list[item]


class QuantNumbersH216O (QuantNumbers):
    def __init__(self, v1, v2, v3, J, Ka, Kc):
        QuantNumbers.__init__(self, [v1, v2, v3, J, Ka, Kc])
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.J = J
        self.Ka = Ka
        self.Kc = Kc


class QuantNumbersN2O (QuantNumbers):
    def __init__(self, v1, v2, J, l):
        QuantNumbers.__init__(self, [v1, v2, J, l])
        self.v1 = v1
        self.v2 = v2
        self.J = J
        self.l = l


class SymType(Enum):
    A1 = 1
    A2 = 2
    B1 = 3
    B2 = 4

    @staticmethod
    def from_int(i):
        return SymType(i)

    @staticmethod
    def from_str(s):
        return SymType(s.upper())


class States:
    def __init__(self, list_states):
        self.__states = list_states
        self.__states = self.__sort()   # always sorted

    def __len__(self):
        return len(self.__states)

    def __iter__(self):
        return iter(self.__states)

    def __getitem__(self, item):
        return self.__states[item]

    def remove_item(self, item):
        self.__states.remove(item)

    def __sort(self):
        return sorted(self.__states, key=lambda x: (x.J, x.sym, x.E))

    def write_to_file(self, out_file_name):
        with open(out_file_name, 'w') as f:
            for state in self.__states:
                if len(state.qn) == 6:
                    f.write("{0:2d} {1:2s} {2:4d} {3:15.6f} {4:3d} {5:2d} {6:2d} {7:5.2f} {8:3d} {9:2d} {10:2d}\n".format(
                        state.J, SymType(state.sym).name, state.N, state.E, state.qn.v1, state.qn.v2, state.qn.v3, state.w,
                        state.qn.J, state.qn.Ka, state.qn.Kc
                    ))
                elif len(state.qn) == 4:
                    f.write("{0:2d} {1:2s} {2:4d} {3:15.6f} {4:3d} {5:2d} {6:2d} {7:5.2f} {8:3d}\n".format(
                        state.J, SymType(state.sym).name, state.N, state.E, state.qn.v1, state.qn.v2, state.qn.J, state.w, state.qn.l
                    ))
                else:
                    print("By now only the formats with 4 and 6 quantum numbers are supported")


class Status(Enum):
    FOUND = 0
    OUTLIER = 1
    NOT_FOUND = 2

    @staticmethod
    def from_int(i):
        return Status(i)

    @staticmethod
    def from_str(s):
        return Status(s.upper())


class ComparedState:
    def __init__(self, E_calc, E_exp, J, sym, N, w, E_diff, status=None, qn=None):
        self.E_calc = E_calc
        self.E_exp = E_exp
        self.E_diff = E_diff
        self.J = J
        self.sym = sym
        self.N = N
        self.w = w
        self.status = status
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
        return sorted(self.__comp_states, key=lambda x: (x.J, x.sym, x.E_exp))

    def write_to_file(self, out_file_name, filter_by=-1.0, mode='w'):
        with open(out_file_name, mode) as f:
            if filter_by < 0.0:
                states_to_write = self.__comp_states
            else:
                # write only the states with w = filter_by
                states_to_write = list(filter(lambda x: abs(x.w - filter_by) < 1e-3, self.__comp_states))

            for state in states_to_write:
                if len(state.qn) == 6:
                    f.write("{0:2d} {1:2s} {2:4d} {3:15.6f} {4:15.6f} {5:15.6f} {6:3d} {7:2d} {8:2d} {9:5.2f} {10:3d} {11:2d} {12:2d} {13:10s}\n".format(
                        state.J, SymType(state.sym).name, state.N, state.E_exp, state.E_calc, state.E_diff,
                        state.qn.v1, state.qn.v2, state.qn.v3, state.w, state.qn.J, state.qn.Ka, state.qn.Kc, Status(state.status).name
                    ))
                elif len(state.qn) == 4:
                    f.write(
                        "{0:2d} {1:2s} {2:4d} {3:15.6f} {4:15.6f} {5:15.6f} {6:3d} {7:2d} {8:2d} {9:5.2f} {10:3d} {11:10s}\n".format(
                            state.J, SymType(state.sym).name, state.N, state.E_exp, state.E_calc, state.E_diff,
                            state.qn.v1, state.qn.v2, state.qn.J, state.w, state.qn.l, Status(state.status).name
                    ))
                else:
                    print("By now only the formats with 4 and 6 quantum numbers are supported")

        return states_to_write

    def parse_file(self, out_file_name):
        comp_states = []
        with open(out_file_name, 'r') as f:
            for line in f:
                list_values = line.strip().split()
                if len(list_values) == 14:
                    qn = QuantNumbersH216O(v1=int(list_values[6]), v2=int(list_values[7]), v3=int(list_values[8]),
                                           J=int(list_values[10]), Ka=int(list_values[11]), Kc=int(list_values[12]))
                    status_txt = list_values[13]
                elif len(list_values) == 12:
                    qn = QuantNumbersN2O(v1=int(list_values[6]), v2=int(list_values[7]), J=int(list_values[8]),
                                         l=int(list_values[10]))
                    status_txt = list_values[11]
                else:
                    print("By now only the formats with 4 and 6 quantum numbers are supported")

                if status_txt == 'FOUND':
                    status = Status.FOUND
                elif status_txt == 'OUTLIER':
                    status = Status.OUTLIER
                else:
                    status = Status.NOT_FOUND

                if list_values[1] == 'A1':
                    sym = SymType.A1
                elif list_values[1] == 'A2':
                    sym = SymType.A2
                elif list_values[1] == 'B1':
                    sym = SymType.B1
                elif list_values[1] == 'B2':
                    sym = SymType.B2
                else:
                    print("Only symmetry types A1, A2, B1, and B2 are supported")

                comp_states.append(ComparedState(float(list_values[4]), float(list_values[3]), int(list_values[0]),
                                                 sym.value, int(list_values[2]), float(list_values[9]), float(list_values[5]),
                                                 status, qn))

        self.__comp_states = comp_states
        return self.__comp_states







