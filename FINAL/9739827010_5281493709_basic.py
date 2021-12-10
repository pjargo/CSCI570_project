import pandas as pd
import numpy as np
import sys
import guppy
from guppy import hpy
import time

input_file = sys.argv[1]
output_file = 'output.txt'

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


# FUNCTIONS

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


# Generate full strings
def generate_string(base, indices):
    cumulative = base
    
    for i in indices:
        cumulative = cumulative[:i + 1] + cumulative + cumulative[i + 1:]
        
    return cumulative
        

# Bottom Up Pass
def alignment(X, Y):
    
    global delta
    global alphas
    
    # size of matrix n x m
    # Y -> has length 1, 2, ...j..., n   ----> Rows
    # X -> has length 1, 2, ...i..., m   ----> Cols
    n = len(Y)    # Rows
    m = len(X)    # Cols

    A = [ [ 0 for i in range(m + 1) ] for j in range(n + 1) ]    # Create empty array
    
    # Initialize A[i, 0] and A[0, j]
    # Matrix access looks like A[row][column] = A[j][i]
    for i in range(m + 1):
        A[0][i] = i * delta

    for j in range(n + 1):
        A[j][0] = j * delta

    # Run recurrence
    for j in range(n + 1)[1:]:
        for i in range(m + 1)[1:]:

            A[j][i] = min(alphas[(X[i - 1], Y[j - 1])] + A[j - 1][i - 1], 
                         delta + A[j][i - 1], 
                         delta + A[j - 1][i])
    
    return (A[n][m], A)


# Top down pass to determine the actual string alignments
def top_down_pass(A, X, Y):
    
    global delta
    global alphas
    
    new_X = ''
    new_Y = ''
    
    j = len(Y)    # Rows
    i = len(X)    # Cols
    
    
    while (i, j) != (0, 0):
        lookback = []    # List to store elements of optimization. We later take minimum from it
        
        # alpha + Opt(i-1, j-1)
        if i == 0 or j == 0:
            lookback.append((float('inf'), 'a_paired'))
        else:
            alpha = alphas[(X[i - 1], Y[j - 1])]
            paired = alpha + A[j - 1][i - 1]
            lookback.append((paired, 'a_paired'))

        # delta + Opt(i, j-1)
        if j == 0:
            lookback.append((float('inf'), 'b_i_gap'))
        else:
            i_gap = delta + A[j - 1][i]
            lookback.append((i_gap, 'b_i_gap'))

        # delta + Opt(i-1, j)
        if i == 0:
            lookback.append((float('inf'), 'c_j_gap'))
        else:
            j_gap = delta + A[j][i - 1]
            lookback.append((j_gap, 'c_j_gap'))

        lookback.sort()

        parent = lookback[0]

        # Update the new strings new_X and new_Y
        if parent[1] == 'a_paired':
            new_X = X[i - 1] + new_X    # this is actuall i but we have to subract one for indexing purposes
            new_Y = Y[j - 1] + new_Y

            i = i - 1    # go to (i - 1, j - 1)
            j = j - 1

        elif parent[1] == 'b_i_gap':
            new_X = '_' + new_X
            new_Y = Y[j - 1] + new_Y

            j = j - 1    # go to (i, j - 1)

        elif parent[1] == 'c_j_gap':
            new_X = X[i - 1] + new_X
            new_Y = '_' + new_Y

            i = i - 1    # go to (i - 1, j)
            
            
    return (new_X, new_Y)
        

# Print Output
def get_output(X, Y, alignment_cost, time_taken, memory):
    global output_file
    
    out_file = open(output_file,"w")
 
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


# Read in file
file_data = read_file(input_file)
base_string_X = file_data['base_string_X']
base_string_Y = file_data['base_string_Y']
X_indices = file_data['X_indices']
Y_indices = file_data['Y_indices']

# Generate Strings
X = generate_string(base_string_X, X_indices)
Y = generate_string(base_string_Y, Y_indices)

# initiate heap
heap = hpy()
heap.setref()

heap_status_before = heap.heap()
#print("Before Heap Size : ", heap_status_before.size * 0.001, " kilobytes\n")

#time
start_time = time.time()

# Do bottom up pass to create alignment matrix
alignment_data = alignment(X, Y)
## A = alignment_data[1]
## alignment_cost = alignment_data[0]

# Top down pass to determine alignment strings
final_strings = top_down_pass(alignment_data[1], X, Y)
## X_new = final_strings[0]
## Y_new = final_strings[1]

end_time = time.time()
time_taken = end_time - start_time

heap_status_after = heap.heap()
# print("After Heap Size : ", heap_status_after.size * 0.001, " kilobytes\n")
memory = (heap_status_after.size - heap_status_before.size) * 0.001


# Print output
get_output(final_strings[0], final_strings[1], alignment_data[0], time_taken, memory)