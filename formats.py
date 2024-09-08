import states


class Format:
    def __init__(self, name):
        self.__name = name

    def parse_file(self):
        pass


class CalcFormat (Format):
    def __init__(self, file_name, Ezero=0.0, au_to_cm=219474.624):
        Format.__init__(self, 'fort.14.fmt')
        self.__file_name = file_name
        self.__Ezero = Ezero
        self.__au_to_cm = au_to_cm

    def parse_file(self):
        list_states = []
        sym = 0
        with open(self.__file_name, 'r') as f:
            line = f.readline()
            while line != '':
                try:
                    J, n_levs = int(line.strip().split()[0]), int(line.strip().split()[5])
                    sym += 1
                    list_ens = []
                    r = n_levs // 4 + 1 if n_levs % 4 != 0 else n_levs // 4
                    for l in range(r):
                        line = f.readline().replace('D', 'E')
                        list_ens.extend(list(map(lambda x: float(x) * self.__au_to_cm - self.__Ezero, line.strip().split())))
                    for en in list_ens:
                        state = states.State(en, J, sym, N=list_ens.index(en) + 1)
                        list_states.append(state)
                except ValueError:
                    print("Something went wrong")

                line = f.readline()

        return states.States(list_states)


class ExpFormat (Format):
    def __init__(self, format_name, file_name, J_list, J_place):
        Format.__init__(self, format_name)
        self.__file_name = file_name
        self.__J_list = J_list
        self.__J_place = J_place
        self.__lines = self.__read_file()

    def __read_file(self):
        with open(self.__file_name, 'r') as f:
            lines = f.readlines()

        read_lines = []
        for line in lines:
            line = line.strip().split()
            read_lines.append(line)

        return read_lines

    def filter_Jlist_lines(self):
        # It's possible to filter out experimental states with J values from the J_list
        return list(filter(lambda x: x[self.__J_place] in self.__J_list, self.__lines))

    def sym_definition(self, qns=None, lbl=''):
        pass


class ExpFormatH216O (ExpFormat):
    def __init__(self, file_name, J_list, J_place):
        ExpFormat.__init__(self, 'exp.h2-16o.fmt', file_name, J_list, J_place)

    def sym_definition(self, qns=None, lbl=''):
        if not qns.J % 2 and not qns.v3 % 2 and not qns.Ka % 2 and not qns.Kc % 2:
            return states.SymType.A1.value
        elif not qns.J % 2 and qns.v3 % 2 and qns.Ka % 2 and not qns.Kc % 2:
            return states.SymType.A1.value
        elif qns.J % 2 and not qns.v3 % 2 and not qns.Ka % 2 and qns.Kc % 2:
            return states.SymType.A1.value
        elif qns.J % 2 and qns.v3 % 2 and qns.Ka % 2 and qns.Kc % 2:
            return states.SymType.A1.value
        elif not qns.J % 2 and not qns.v3 % 2 and not qns.Ka % 2 and qns.Kc % 2:
            return states.SymType.A2.value
        elif not qns.J % 2 and qns.v3 % 2 and qns.Ka % 2 and qns.Kc % 2:
            return states.SymType.A2.value
        elif qns.J % 2 and not qns.v3 % 2 and not qns.Ka % 2 and not qns.Kc % 2:
            return states.SymType.A2.value
        elif qns.J % 2 and qns.v3 % 2 and qns.Ka % 2 and not qns.Kc % 2:
            return states.SymType.A2.value
        elif not qns.J % 2 and not qns.v3 % 2 and qns.Ka % 2 and not qns.Kc % 2:
            return states.SymType.B1.value
        elif not qns.J % 2 and qns.v3 % 2 and not qns.Ka % 2 and not qns.Kc % 2:
            return states.SymType.B1.value
        elif qns.J % 2 and not qns.v3 % 2 and qns.Ka % 2 and qns.Kc % 2:
            return states.SymType.B1.value
        elif qns.J % 2 and qns.v3 % 2 and not qns.Ka % 2 and qns.Kc % 2:
            return states.SymType.B1.value
        elif not qns.J % 2 and not qns.v3 % 2 and qns.Ka % 2 and qns.Kc % 2:
            return states.SymType.B2.value
        elif not qns.J % 2 and qns.v3 % 2 and not qns.Ka % 2 and qns.Kc % 2:
            return states.SymType.B2.value
        elif qns.J % 2 and not qns.v3 % 2 and qns.Ka % 2 and not qns.Kc % 2:
            return states.SymType.B2.value
        elif qns.J % 2 and qns.v3 % 2 and not qns.Ka % 2 and not qns.Kc % 2:
            return states.SymType.B2.value

    def parse_file(self):
        filtered_lines = self.filter_Jlist_lines()
        list_states = []
        for line in filtered_lines:
            # 0 ~ v1  1 ~ v2  2 ~ v3  3 ~ J  4 ~ Ka  5 ~ Kc  6 ~ E
            qns = states.QuantNumbersH216O(v1=int(line[0]), v2=int(line[1]), v3=int(line[2]), J=int(line[3]), Ka=int(line[4]), Kc=int(line[5]))
            sym = self.sym_definition(qns)
            state = states.State(float(line[6]), int(line[3]), sym, qn=qns)
            list_states.append(state)

        return states.States(list_states)


class ExpFormatN2O (ExpFormat):
    def __init__(self, file_name, J_list, J_place):
        ExpFormat.__init__(self, 'exp.n2o.fmt', file_name, J_list, J_place)

    def sym_definition(self, qns=None, lbl=''):
        if lbl == 'e':
            return states.SymType.A1.value
        elif lbl == 'f':
            return states.SymType.A2.value
        else:
            print("Strange symmetry label got!")

    def parse_file(self):
        filtered_lines = self.filter_Jlist_lines()
        list_states = []
        for line in filtered_lines:
            # 0 ~ q1  1 ~ q2  2 ~ q3  3 ~ sym  4 ~ J  5 ~ E
            qns = states.QuantNumbersN2O(v1=int(line[0]), v2=int(line[2]), l=int(line[1]), J=int(line[4]))
            sym = self.sym_definition(lbl=line[3])
            state = states.State(float(line[5]), int(line[4]), sym, qn=qns)
            list_states.append(state)

        return states.States(list_states)



