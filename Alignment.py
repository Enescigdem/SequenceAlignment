import numpy as np

import argparse


def matrix(m, s1, s2):
    c = 0
    print("", end="\t\t")
    for l in s2:
        print(l, end="\t")
    print()
    for i in range(m.shape[0]):
        if not c == 0:
            print(s1[i - 1], end="\t")
        else:
            print("", end="\t")
        for j in range(m.shape[1]):
            print(str(m[i][j]), end="\t")
        print()

        c += 1


def create_scorematrixdict(scorematrixpath):
    scorefile = open(scorematrixpath, 'r')
    scorematrix_dict = {}
    columns = []
    counter = 0
    for line in scorefile.readlines():
        if line.startswith("#"):
            continue
        row = line.split()
        if counter == 0:
            columns = row
            for elt in row:
                eltdict = {}
                scorematrix_dict[elt] = eltdict
        else:

            for column in columns:
                scorematrix_dict[row[0]][column] = row[columns.index(column) + 1]

        counter += 1
    scorefile.close()
    return scorematrix_dict


def get_scoring_score(sdict, seq1_element, seq2_element):
    if seq1_element in sdict.keys():
        if seq2_element in sdict.keys():
            return sdict[seq1_element][seq2_element]

        return sdict[seq1_element]["*"]
    else:
        if seq2_element in sdict.keys():
            return sdict[seq2_element]["*"]
        else:
            return sdict["*"]["*"]


def get_max_score(matrix, ii, jj, scoredict, opening, extension, seq1, seq2, method):
    diagonal = matrix[ii - 1][jj - 1] + np.int16(get_scoring_score(scoredict, seq1[ii - 1], seq2[jj - 1]))
    fromleft = matrix[ii][jj - 1] + opening
    fromtop = matrix[ii - 1][jj] + extension

    scores = (diagonal, fromtop, fromleft)
    if method == "global":
        return max(scores), np.argmax(scores)
    elif method == "local":

        return max(0, diagonal, fromleft, fromtop), np.argmax(scores)


def needlemanwunsch(sequence1, sequence2, scoredict, gap_open_penalty, gap_extension_penalty,outfile):
    alignment_matrix = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
    row_c, col_c = alignment_matrix.shape
    outfile.write("Needleman-Wunsch\n")
    trace = np.full((len(sequence1) + 1, len(sequence2) + 1), 0)
    for i in range(len(sequence2) + 1):
        trace[0, i] = 2
    for i in range(len(sequence1) + 1):
        trace[i, 0] = 1

    for row in range(row_c):
        alignment_matrix[row][0] = gap_open_penalty * (row)
    for col in range(col_c):
        alignment_matrix[0][col] = gap_extension_penalty * (col)
    for i in range(1, row_c):
        for j in range(1, col_c):
            alignment_matrix[i][j], tracedirection = get_max_score(alignment_matrix, i, j, scoredict, gap_open_penalty,
                                                                   gap_extension_penalty, sequence1, sequence2,
                                                                   "global")
            trace[i, j] = tracedirection

    i = len(sequence1)
    j = len(sequence2)
    firstline = []
    match_mismatchline = []
    thirdline = []
    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == 0:
            firstline.append(sequence1[i - 1])
            thirdline.append(sequence2[j - 1])
            i -= 1
            j -= 1
        if direction == 1:
            firstline.append(sequence1[i - 1])
            thirdline.append("-")
            i -= 1
        if direction == 2:
            firstline.append("-")
            thirdline.append(sequence2[j - 1])
            j -= 1
    match_count = 0
    summary = ""
    for i, j in zip(firstline, thirdline):
        if i != j:
            match_mismatchline.append(" ")
            summary += "N"
        else:
            match_mismatchline.append("|")
            summary += "M"
        if i == j:
            match_count += 1
    indel_num = 0
    match_num = 0
    for i in range(len(summary)):

        if summary[i] == "M":
            match_num += 1
        elif summary[i] == "N":
            indel_num += 1
    firstline = list(reversed(firstline))
    match_mismatchline = list(reversed(match_mismatchline))
    thirdline = list(reversed(thirdline))
    outfile.write("".join(firstline))
    outfile.write("\n")
    outfile.write("".join(match_mismatchline))
    outfile.write("\n")
    outfile.write("".join(thirdline))
    outfile.write("\n")
    print("".join(firstline))
    print("".join(match_mismatchline))
    print("".join(thirdline))
    i, j = len(sequence1), len(sequence2)
    outfile.write("Raw alignment score  :{} \n".format(alignment_matrix[-1][-1]))
    print("Raw alignment score  :{} \n".format(alignment_matrix[-1][-1]))
    outfile.write("Match rate :{} percent \n".format(match_count * 100 / len(firstline)))
    print("Match rate :{} percent".format(match_count * 100 / len(firstline)))
    outfile.write("{} matches {} insertion-deletions in aligned sequence \n".format(match_num, indel_num))
    print("{} matches {} insertion-deletions in aligned sequence".format(match_num, indel_num))





def smithwaterman(sequence1, sequence2, scoredict, gap_open_penalty, gap_extension_penalty,outfile):
    alignment_matrix = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
    row_c, col_c = alignment_matrix.shape
    trace = np.full((len(sequence1) + 1, len(sequence2) + 1), 0)
    outfile.write("Smith-Waterman\n")
    for row in range(row_c):
        alignment_matrix[row][0] = 0
    for col in range(col_c):
        alignment_matrix[0][col] = 0

    for i in range(1, row_c):
        for j in range(1, col_c):
            alignment_matrix[i][j], trace_index = get_max_score(alignment_matrix, i, j, scoredict, gap_open_penalty,
                                                                gap_extension_penalty, sequence1, sequence2, "local")

            trace[i, j] = trace_index
    max = 0
    for irow in range(len(alignment_matrix)):
        for icol in range(len(alignment_matrix[0])):
            if alignment_matrix[irow][icol] > max:
                max = alignment_matrix[irow][icol]
                i, j = irow, icol

    firstline = []
    match_mismatchline = []
    thirdline = []
    while i > 0 or j > 0:
        direction = trace[i][j]
        if alignment_matrix[i, j] == 0:
            break
        if direction == 0:  # diagonal
            firstline.append(sequence1[i - 1])
            thirdline.append(sequence2[j - 1])

            i -= 1
            j -= 1
        elif direction == 1:  # from left
            firstline.append(sequence1[i - 1])
            thirdline.append("-")

            i -= 1
        elif direction == 2:  # from top
            firstline.append("-")
            thirdline.append(sequence2[j - 1])
            j -= 1
    local_match_count = 0
    summary = ""
    for i, j in zip(firstline, thirdline):
        if i == j:
            local_match_count += 1
        if i != j:
            match_mismatchline.append(" ")
            summary += "N"
        else:
            match_mismatchline.append("|")
            summary += "M"
    indel_num = 0
    match_num = 0
    for i in range(len(summary)):

        if summary[i] == "M":
            match_num += 1
        elif summary[i] == "N":
            indel_num += 1

    firstline = list(reversed(firstline))
    match_mismatchline = list(reversed(match_mismatchline))
    thirdline = list(reversed(thirdline))
    print("".join(firstline))
    print("".join(match_mismatchline))
    print("".join(thirdline))
    outfile.write("".join(firstline))
    outfile.write("\n")
    outfile.write("".join(match_mismatchline))
    outfile.write("\n")
    outfile.write("".join(thirdline))
    outfile.write("\n")
    outfile.write("Raw alignment score  :{} \n".format(max))
    print("Raw alignment score  :{} \n".format(max))
    outfile.write("Match rate :{} percent \n".format(local_match_count * 100 / len(firstline)))
    print("Match rate :{} percent".format(local_match_count * 100 / len(firstline)))
    outfile.write("{} matches {} insertion-deletions in aligned sequence \n".format(match_num, indel_num))
    print("{} matches {} insertion-deletions in aligned sequence".format(match_num, indel_num))



def main():


    parser = argparse.ArgumentParser()
    parser.add_argument('--method', type=str, default="local", help='global or local')
    parser.add_argument('--scoringmatrix_filepath', type=str, default="IDENTITY.txt", help='Matrix txt path to be used')
    parser.add_argument('--input_path', type=str, default="alignment_input.txt",
                        help='txt file path of input sequences  to be used')

    parser.add_argument('--gap_extension_penalty', type=int, default=-8)
    parser.add_argument('--gap_open_penalty', type=int, default=-8)
    parser.add_argument('--outfile', type=str, default="Different_parameter_Alignments.txt")
    opt = parser.parse_args()
    output = open(opt.outfile, "a")
    output.write(opt.method+"\n")
    output.write(opt.scoringmatrix_filepath+"\n")

    output.write("{} gap_open_penalty |  {} gap_extension_penalty\n".format(opt.gap_open_penalty, opt.gap_extension_penalty))
    output.write("\n")
    print(opt)
    inputfile = open(opt.input_path, 'r')
    scoredict = create_scorematrixdict(opt.scoringmatrix_filepath)
    line_count = len(inputfile.readlines())
    inputfile.seek(0)
    for c in range(0, line_count, 2):
        s1 = inputfile.readline()
        s2 = inputfile.readline()
        if opt.method == 'local':
            smithwaterman(s1, s2, scoredict, opt.gap_open_penalty, opt.gap_extension_penalty,output)
        elif opt.method == 'global':
            needlemanwunsch(s1, s2, scoredict, opt.gap_open_penalty, opt.gap_extension_penalty,output)

    output.write("--------------------------------------------------------------------------------------\n")
    output.close()

main()
