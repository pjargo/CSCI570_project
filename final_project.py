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

    A = [ [ 0 for i in range(m) ] for j in range(n) ]    # Create empty array
    
    # Initialize A[i, 0] and A[0, j]
    # Matrix access looks like A[row][column] = A[j][i]
    for i in range(m):
        A[0][i] = i * delta

    for j in range(n):
        A[j][0] = j * delta
    
    # Run recurrence
    for j in range(n)[1:]:
        for i in range(m)[1:]:

            A[j][i] = min(alphas[(X[i], Y[j])] + A[j - 1][i - 1], 
                         delta + A[j][i - 1], 
                         delta + A[j - 1][i])
    
    return A[n - 1][m - 1]


if __name__ == '__main__':
    
    input_file = sys.argv[1]
    output_file = 'output.txt'
    
    file_data = read_file(input_file)
    base_string_X = file_data['base_string_X']
    base_string_Y = file_data['base_string_Y']
    X_indices = file_data['X_indices']
    Y_indices = file_data['Y_indices']
    
    
    
    s1 = "ACTG"     # "ACACTGACTACTGACTGGTGACTACTGACTGG"
    j = [3,6,1]
    s2 = "TACG"     # "TATTATACGCTATTATACGCGACGCGGACGCG"
    k = [1, 2, 9]
    new_s1 = string_generator(string=s1, jk=j)

    print(new_s1)

