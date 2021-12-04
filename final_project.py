"""
Project Name: CSCI 570 Final Project
Author: Luke Nelson, Peter Argo
"""

import sys


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
    # size of matrix n x m
    # Y -> has length 1, 2, ...j..., n   ----> Rows
    # X -> has length 1, 2, ...i..., m   ----> Cols
    n = len(Y)    # Rows
    m = len(X)    # Cols

    A = [ [ 0 for i in range(m) ] for j in range(n) ]    # Create empty array
    
    # Initialize A[i, 0] and A[0, j]
    # Matrix access looks like A[row][column] = A[j][i]
    for i in range(m):
        A[0][i] = i * delta

    for j in range(n):
        A[j][0] = j * delta
    
    # Run recurrence
    for j in range(n):
        for i in range(m):

            A[j][i] = min(alphas[(X[i], Y[j])] + A[j - 1][i - 1], 
                         delta + A[j][i - 1], 
                         delta + A[j - 1][i])
    
    return A[n - 1][m - 1]


if __name__ == '__main__':
    
    input_file = sys.argv[1]
    output_file = 'output.txt'
    
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
    
    
    
    s1 = "ACTG"     # "ACACTGACTACTGACTGGTGACTACTGACTGG"
    j = [3,6,1]
    s2 = "TACG"     # "TATTATACGCTATTATACGCGACGCGGACGCG"
    k = [1, 2, 9]
    new_s1 = string_generator(string=s1, jk=j)

    print(new_s1)
