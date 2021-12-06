"""
Project Name: CSCI 570 Final Project
Author: Luke Nelson, Peter Argo
"""

import sys

delta = 30
alphas = {('A', 'A'): 0,
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


def read_file(input_file):
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

    return {'base_string_X': base_string_X,
           'base_string_Y': base_string_Y,
           'X_indices': X_indices,
           'Y_indices': Y_indices}


def string_generator(string: str, jk: list):
    """
    Gererate the stings from the input files
    :param string: input string
    :param jk: list of multipliers
    :return:
    """
    if len(jk) == 1:
        return string[0:jk[0]+1] + string + string[jk[0]+1::]

    new_string = string[0:jk[0]+1] + string + string[jk[0]+1::]
    return string_generator(new_string, jk[1::])


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


def alignment(X, Y):

    global delta
    global alphas

    # size of matrix n x m
    # Y -> has length 1, 2, ...j..., n   ----> Rows
    # X -> has length 1, 2, ...i..., m   ----> Cols
    n = len(Y)    # Rows
    m = len(X)    # Cols
    print("length of Y string: (n) ", n)
    print("length of X string: (m) ", m)
    A = [ [ 0 for i in range(m + 1) ] for j in range(n + 1) ]    # Create empty array
    T = [[0 for i in range(m + 1)] for j in range(n + 1)]
    # Initialize A[i, 0] and A[0, j]
    # Matrix access looks like A[row][column] = A[j][i]
    for i in range(m + 1):
        A[0][i] = i * delta

    for j in range(n + 1):
        A[j][0] = j * delta

    # Run recurrence
    for j in range(n + 1)[1:]:
        for i in range(m + 1)[1:]:

            diag = alphas[(X[i - 1], Y[j - 1])] + A[j - 1][i - 1]
            left = delta + A[j][i - 1]
            up = delta + A[j - 1][i]

            A[j][i] = min(alphas[(X[i - 1], Y[j - 1])] + A[j - 1][i - 1],
                         delta + A[j][i - 1],
                         delta + A[j - 1][i])

            best_choice = min(diag, left, up)

            if best_choice == diag:
                # T[j][i] = str(i-1) + " " + str(j-1)
                T[j][i] = str(j - 1) + " " + str(i - 1)
            elif best_choice == left:
                # T[j][i] = str(i - 1) + " " + str(j)
                T[j][i] = str(j - 1) + " " + str(i)
            elif best_choice == up:
                # T[j][i] = str(i) + " " + str(j-1)
                T[j][i] = str(j) + " " + str(i - 1)

            new_T = [row[1::] for row in T][1::]

    return (new_T, A)


def reconstruct_soln(matrix, m, n, xstring, ystring):
    """

    :param matrix: The alignment matrix
    :param m:
    :param n:
    :param xstring:
    :param ystring:
    :return:
    """
    global solution
    global solution_stringx
    global solution_stringy
    global x_prev
    global y_prev
    global delta
    global alphas

    if m == 0 and n == 0:
        return

    left = matrix[n][m-1]
    up = matrix[n-1][m]
    diag = matrix[n-1][m-1]

    if matrix[n][m] == diag + alphas[(xstring[m-1], ystring[n-1])]:
        solution_stringx = xstring[m-1] + solution_stringx
        solution_stringy = ystring[n - 1] + solution_stringy
        return reconstruct_soln(matrix=matrix, m=m-1, n=n-1, xstring=xstring, ystring=ystring)

    elif matrix[n][m] == left + delta:
        solution_stringx = xstring[m-1] + solution_stringx
        solution_stringy = "_" + solution_stringy
        return reconstruct_soln(matrix, m-1, n, xstring, ystring)

    elif matrix[n][m] == up + delta:
        solution_stringx = "_" + solution_stringx
        solution_stringy = ystring[n - 1] + solution_stringy
        return reconstruct_soln(matrix, m, n-1, xstring, ystring)


if __name__ == '__main__':
    # input_file = sys.argv[1]
    input_file = "input1.txt"
    output_file = 'output.txt'
    solution = list()
    solution_stringx = ""
    solution_stringy = ""

    file_data = read_file(input_file)
    base_string_X = file_data['base_string_X']
    base_string_Y = file_data['base_string_Y']
    X_indices = file_data['X_indices']
    Y_indices = file_data['Y_indices']

    s1 = "ACTG"
    j = [3,6,1,1]
    s2 = "TACG"
    k = [1, 2, 9, 2]
    # x_prev = len(s1)
    # y_prev = len(s2)

    new_s1 = string_generator(string=s1, jk=j)
    new_s2 = string_generator(string=s2, jk=k)
    print("string 1: ", new_s1)
    print("string 2: ", new_s2)

    x_prev = len(new_s1)
    y_prev = len(new_s2)

    mm = alignment(new_s1, new_s2)
    # mm = alignment(s1, s2)

    reconstruct_soln(matrix=mm[1], m=len(new_s1), n=len(new_s2), xstring=new_s1, ystring=new_s2)
    # reconstruct_soln(matrix=mm[1], m=len(s1), n=len(s2), xstring=s1, ystring=s2)

    print(solution[::-1])
    print("our solution for string 1(x): ", solution_stringx[0:50], solution_stringx[len(solution_stringx)-50::])
    print("our solution for string 1(y): ", solution_stringy[0:50], solution_stringy[len(solution_stringy)-50::])

    print("output1.txt for string 1(x):  ", "_A_CA_CACT__G__A_C_TAC_TGACTG_GTGA__C_TACTGACTGGAC", "GTGA__C_TACTGACTGGACTGACTACTGACTGGTGACTACT_GACTG_G")
    print("output1.txt for string 1(y):  ", "TATTATTA_TACGCTATTATACGCGAC_GCG_GACGCGTA_T_AC__G_C", "G_GACGCGTA_T_AC__G_CT_ATTA_T_AC__GCGAC_GC_GGAC_GCG")
