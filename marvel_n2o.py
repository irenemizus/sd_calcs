import sys


def parse_exp_file(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()

    parsed_lines = []
    for line in lines:
        line = line.strip().split()
        parsed_lines.append(line)

    return parsed_lines


def filter_Jlist_lines(lines, J_list, J_place):
    return list(filter(lambda x: x[J_place] in J_list, lines))


def sym_definition(sym):
    if sym == 'e':
        return 'A1'
    elif sym == 'f':
        return 'A2'


def get_eners(file_name, J_list):
    lines = parse_exp_file(file_name)
    filtered_lines = filter_Jlist_lines(lines, J_list, 4)
    new_lines = []
    for line in filtered_lines:
        # 0 ~ q1  1 ~ q2  2 ~ q3  3 ~ sym  4 ~ J  5 ~ E
        new_line = [int(line[4]), sym_definition(line[3]), float(line[5])]
        if new_line[1] == 'A1':
            new_line[1] = '1'
        elif new_line[1] == 'A2':
            new_line[1] = '2'
        elif new_line[1] == 'B1':
            new_line[1] = '3'
        elif new_line[1] == 'B2':
            new_line[1] = '4'

        new_lines.append(new_line)

    return sorted(new_lines, key=lambda x: (x[0], x[1], x[2]))


if __name__ == '__main__':
    lines = parse_exp_file(sys.argv[1])
    J_list = ['2', '15', '20']
    filtered_lines = filter_Jlist_lines(lines, J_list, 4)

    with open(file=sys.argv[2], mode='w') as f:
        new_lines = []
        for line in filtered_lines:
            # 0 ~ q1  1 ~ q2  2 ~ q3  3 ~ sym  4 ~ J  5 ~ E
            new_line = [int(line[4]), sym_definition(line[3]), '0', float(line[5]), line[0], line[1], line[2], ' 1.00']
            new_lines.append(new_line)

        sorted_list = sorted(new_lines, key=lambda x: (x[0], x[1], x[3]))
        for J in J_list:
            for sym in ['A1', 'A2']:
                N = 0
                for line in sorted_list:
                    if line[0] == int(J) and line[1] == sym:
                        N += 1
                        line[2] = str(N)

        for line in sorted_list:
            #f.write('\t'.join(line) + '\n')
            f.write("{0:2d} {1:2s} {2:4d} {3:15.6f} {4:3d} {5:2d} {6:2d} {7:5s}\n".format(
                int(line[0]), line[1], int(line[2]), line[3], int(line[4]), int(line[5]), int(line[6]), line[7]
            ))
            print(line)

        print(len(filtered_lines), len(sorted_list))



