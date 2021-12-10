"""
Project Name: CSCI 570 Final Project
Author: Luke Nelson, Peter Argo
"""

import sys
import time
from guppy import hpy

def generate_string(string: str, jk: list):
    """
    Gererate the stings from the input files
    :param string: input string
    :param jk: list of multipliers
    :return:
    """
    if len(jk) == 1:
        return string[0:jk[0] + 1] + string + string[jk[0] + 1::]

    new_string = string[0:jk[0] + 1] + string + string[jk[0] + 1::]
    return generate_string(new_string, jk[1::])


class sequenceAlign():
    def __init__(self):
        self.delta = 30
        self.input_data = dict()
        self.alphas = {('A', 'A'): 0,
                      ('A', 'C'): 110,
                      ('A', 'G'): 48,
                      ('A', 'T'): 94,
                      ('C', 'A'): 110,
                      ('C', 'C'): 0,
                      ('C', 'G'): 118,
                      ('C', 'T'): 48,
                      ('G', 'A'): 48,
                      ('G', 'C'): 118,
                      ('G', 'G'): 0,
                      ('G', 'T'): 110,
                      ('T', 'A'): 94,
                      ('T', 'C'): 48,
                      ('T', 'G'): 110,
                      ('T', 'T'): 0}

        self.solution_stringx = ""
        self.solution_stringy = ""
        self.total_score_dnc = 0

    def read_file(self, input_file):
        # Open input file and determine base_string_X with its indices and base_string_Y with its indices
        txt = open(input_file, "r")

        first_line = True
        X_switch = True

        base_string_X = None
        base_string_Y = None
        X_indices = []
        Y_indices = []

        for line in txt:
            line = line.rstrip()

            if first_line == True:
                base_string_X = line
                first_line = False

            elif X_switch == True:
                try:
                    line = int(line)
                    X_indices.append(line)
                except:
                    base_string_Y = line
                    X_switch = False
            else:
                try:
                    Y_indices.append(int(line))
                except:
                    break

        txt.close()

        self.input_data = {'base_string_X': base_string_X,
                             'base_string_Y': base_string_Y,
                             'X_indices': X_indices,
                             'Y_indices': Y_indices}

        return {'base_string_X': base_string_X,
                'base_string_Y': base_string_Y,
                'X_indices': X_indices,
                'Y_indices': Y_indices}

    def alignment(self, X, Y):
        """
        Primary alignmen function using dynamic programing
        :param X: input string along x-axis
        :param Y: input string along y-axis
        :return: optimal alignment matrix of the two sctrings
        """

        # size of matrix n x m
        # Y -> has length 1, 2, ...j..., n   ----> Rows
        # X -> has length 1, 2, ...i..., m   ----> Cols
        n = len(Y)  # Rows
        m = len(X)  # Cols

        A = [[0 for i in range(m + 1)] for j in range(n + 1)]  # Create empty array
        # Initialize A[i, 0] and A[0, j]
        # Matrix access looks like A[row][column] = A[j][i]
        for i in range(m + 1):
            A[0][i] = i * self.delta

        for j in range(n + 1):
            A[j][0] = j * self.delta

        # Run recurrence
        for j in range(n + 1)[1:]:
            for i in range(m + 1)[1:]:
                A[j][i] = min(self.alphas[(X[i - 1], Y[j - 1])] + A[j - 1][i - 1],
                              self.delta + A[j][i - 1],
                              self.delta + A[j - 1][i])

        return A

    def reconstruct_soln(self, matrix, m, n, xstring, ystring):
        """

        :param matrix: The alignment matrix
        :param m: length of xstring
        :param n: length of y string
        :param xstring:
        :param ystring:
        :return: string of letters in the xstring and ystring determined from the matrix
        """

        if m == 0 or n == 0:
            if m == 0 and n == 0:  # Done
                tempx, tempy = self.solution_stringx, self.solution_stringy
                self.reset_solution_stings()
                return tempx, tempy
            elif n == 0:  # Go left
                if matrix[n][m] != matrix[n][m - 1] + self.delta:
                    print("Something needs to be debugged")
                self.solution_stringx = xstring[m - 1] + self.solution_stringx
                self.solution_stringy = "_" + self.solution_stringy
                return self.reconstruct_soln(matrix=matrix, m=m - 1, n=n, xstring=xstring, ystring=ystring)
            else:  # Go Up
                if matrix[n][m] != matrix[n - 1][m] + self.delta:
                    print("Something needs to be debugged")
                self.solution_stringx = "_" + self.solution_stringx
                self.solution_stringy = ystring[n - 1] + self.solution_stringy
                return self.reconstruct_soln(matrix=matrix, m=m, n=n - 1, xstring=xstring, ystring=ystring)

        if m != 0 and n != 0:

            left = matrix[n][m - 1]
            up = matrix[n - 1][m]
            diag = matrix[n - 1][m - 1]
            # print(diag + self.alphas[(xstring[m-1], ystring[n-1])], left + self.delta, up + self.delta)

            if matrix[n][m] == diag + self.alphas[(xstring[m - 1], ystring[n - 1])]:
                self.solution_stringx = xstring[m - 1] + self.solution_stringx
                self.solution_stringy = ystring[n - 1] + self.solution_stringy
                return self.reconstruct_soln(matrix=matrix, m=m - 1, n=n - 1, xstring=xstring, ystring=ystring)

            if matrix[n][m] == up + self.delta:
                self.solution_stringx = "_" + self.solution_stringx
                self.solution_stringy = ystring[n - 1] + self.solution_stringy
                return self.reconstruct_soln(matrix=matrix, m=m, n=n - 1, xstring=xstring, ystring=ystring)

            if matrix[n][m] == left + self.delta:
                self.solution_stringx = xstring[m - 1] + self.solution_stringx
                self.solution_stringy = "_" + self.solution_stringy
                return self.reconstruct_soln(matrix=matrix, m=m - 1, n=n, xstring=xstring, ystring=ystring)

    def dnc_alignment(self, xstring, ystring, a, b, c, d):
        """
        Implement the divide and conquer memory efficient solution to the sequence alignment problem

        :param xstring: input string along x-axis
        :param ystring: input string along y-axis
        :param a: start index for xstring
        :param b: end index for xstring
        :param c: start index for ystring
        :param d: end index for ystring
        :return:
        """
        m = b - a
        n = d - c

        # Base Case
        if m <= 1 or n <= 1:  # Confirmed using or not and
            M = self.alignment(xstring[a:b], ystring[c:d])  # Get the similarity matrix of the small matrix
            score = M[len(M) - 1][len(M[0]) - 1]    # Get the overall score of the submatrix
            string_finalx, string_finaly = self.reconstruct_soln(matrix=M, m=m, n=n, xstring=xstring[a:b],
                                                                 ystring=ystring[c:d])
            return string_finalx, string_finaly, score

        # Divide: Split xstring into two substrings
        position_i = (a + b) // 2

        # Conquer:
        # Determine where to split the y string - find the optimal point
        # compute matrices for left string and right string and find the best score between the two sub-problems

        a1 = self.get_optimal_pt(xstring[a:position_i], ystring[c:d])
        a2 = self.get_optimal_pt(xstring[position_i:b][::-1], ystring[c:d][::-1])[::-1]  # Reverse the right side

        # Divide:
        # Get the optimal index from vectors a1 and a1. This point is where we divide the ystring
        position_j = 0
        min_val = float("inf")
        for j in range(1, len(a1)):  # range(len(a1))??
            score = a1[j][0] + a2[j][0]
            if score < min_val:
                min_val = score
                position_j = c + j

        # Recursive call to the function for each subproblem
        s1x, s1y, score1 = self.dnc_alignment(xstring=xstring, ystring=ystring, a=a, b=position_i, c=c, d=position_j)
        s2x, s2y, score2 = self.dnc_alignment(xstring=xstring, ystring=ystring, a=position_i, b=b, c=position_j, d=d)

        self.set_total_score_dnc(score=score1 + score2)

        return s1x + s2x, s1y + s2y, score1 + score2

    def get_optimal_pt(self, seqx, seqy):
        """
        Get the last column score values which are needed to get the optimal point to divide the ystring. The last
        Dynnamic programming to get the score values
        After computing the score for a column, update the previous column with values we just calculated. This
        saves on memory because we only need to get the last column values

        :param seqx:
        :param seqy:
        :return: a 2 x n vecotor with the optimal scores to align seqx and seqy
        """

        m = len(seqx)
        n = len(seqy)
        a = list()

        # Initialize our data structure
        for j in range(n + 1):
            a.append([j * self.delta, 0])

        # Set the zero and delta values, then compute the scores
        for i in range(1, m + 1):
            a[0][1] = i * self.delta

            for j in range(1, n + 1):
                a[j][1] = min(self.alphas[(seqx[i - 1], seqy[j - 1])] + a[j - 1][0],
                              self.delta + a[j][0],
                              self.delta + a[j - 1][1])

            # Shift column 2 values into column 1 values
            for j in range(n + 1):
                a[j][0] = a[j][1]

        return a

    def reset_solution_stings(self):
        self.solution_stringx = ""
        self.solution_stringy = ""

    def set_total_score_dnc(self, score):
        self.total_score_dnc = score

    def get_total_score_dnc(self):
        return self.total_score_dnc

    @staticmethod
    def get_output(output_file, X, Y, alignment_cost, time_taken, memory):
        out_file = open(output_file, "w")

        out_file.write(X[:50] + ' ' + X[-50:])
        out_file.write('\n')
        out_file.write(Y[:50] + ' ' + Y[-50:])
        out_file.write('\n')
        out_file.write(str(float(alignment_cost)))
        out_file.write('\n')
        out_file.write(str(time_taken))
        out_file.write('\n')
        out_file.write(str(memory))
        out_file.write('\n')

        out_file.close()

    @staticmethod
    def str_length_validation(base_string, full_string, indices):
        full_str_length = len(full_string)
        base_str_length = len(base_string)
        j = len(indices)

        validation_length = (2 ** j) * base_str_length

        if full_str_length == validation_length:
            print('TRUE')
            print('Full String Length: ', full_str_length)
            print('Validation Length: ', validation_length)

        else:
            print('FALSE')
            print('Full String Length: ', full_str_length)
            print('Validation Length: ', validation_length)


if __name__ == '__main__':
    input_file = sys.argv[1]
    # input_file = "input2.txt"
    output_file = 'output.txt'

    sequence_align = sequenceAlign()

    file_data = sequence_align.read_file(input_file)
    base_string_X = file_data['base_string_X']
    base_string_Y = file_data['base_string_Y']
    X_indices = file_data['X_indices']
    Y_indices = file_data['Y_indices']

    new_s1 = generate_string(string=base_string_X, jk=X_indices)
    new_s2 = generate_string(string=base_string_Y, jk=Y_indices)
    # print("string 1: ", new_s1)
    # print("string 2: ", new_s2)

    # alignment_matrix = sequence_align.alignment(X=new_s1, Y=new_s2)
    # total_score1 = alignment_matrix[len(alignment_matrix) - 1][len(alignment_matrix[0]) - 1]
    # print(f"total alignment score for {input_file} data: ", alignment_matrix[len(alignment_matrix)-1][len(alignment_matrix[0])-1])

    # solution_stringx, solution_stringy = sequence_align.reconstruct_soln(matrix=alignment_matrix, m=len(new_s1),
    #                                                                      n=len(new_s2), xstring=new_s1, ystring=new_s2)

    # print()

    # print("simple solution for string 1(x): ", solution_stringx[0:50], solution_stringx[len(solution_stringx)-50::])
    # print("our solution for string 1(y): ", solution_stringy[0:50], solution_stringy[len(solution_stringy)-50::])
    # print()
    # print("output1.txt for string 1(x):  ", "_A_CA_CACT__G__A_C_TAC_TGACTG_GTGA__C_TACTGACTGGAC", "GTGA__C_TACTGACTGGACTGACTACTGACTGGTGACTACT_GACTG_G")
    # print("output1.txt for string 1(y):  ", "TATTATTA_TACGCTATTATACGCGAC_GCG_GACGCGTA_T_AC__G_C", "G_GACGCGTA_T_AC__G_CT_ATTA_T_AC__GCGAC_GC_GGAC_GCG")
    # print("output2.txt for string 1(x):  ", "ACACACTGACTACTGACTGGTGACTACTGACTGGACTGACTACTGACTGG", "CTGGTGACTACTGACTGGACTGACTACTGACTGGTGACTAC_TGACTG_G")
    # print("output2.txt for string 1(y):  ", "__________T__T_A_T__T_A_TAC_G_C__GAC_G____C_GA_T__", "_T__T_A_TAC_G_C__GAC_G____C_GA_T__T_A_TACGCGAC_GCG")


    # =========================== DIVIDE AND CONQUER =======================================
    # Implement the Memory efficient solution Divide and Conquer
    # initiate heap
    heap = hpy()
    heap.setref()
    heap_status_before = heap.heap()
    # start time
    start_time = time.time()

    sequence_align.reset_solution_stings()  # Need to call this each time we run the algo after the first time
    final_stringx, final_stringy, total_score = sequence_align.dnc_alignment(xstring=new_s1, ystring=new_s2,
                                                                             a=0, b=len(new_s1), c=0, d=len(new_s2))

    end_time = time.time()
    time_taken = end_time - start_time
    heap_status_after = heap.heap()
    memory = (heap_status_after.size - heap_status_before.size) * 0.001

    sequence_align.get_output(output_file=output_file, X=final_stringx, Y=final_stringx,
                              alignment_cost=total_score, time_taken=time_taken, memory=memory)

    # print()
    # print("our solution for string 1(x): ", final_stringx[0:50], final_stringx[len(final_stringx)-50::])
    # print("our solution for string 1(y): ", final_stringy[0:50], final_stringy[len(final_stringx)-50::])
    # print(f"total alignment score for {input_file} data using dnc: ", sequence_align.get_total_score_dnc())
