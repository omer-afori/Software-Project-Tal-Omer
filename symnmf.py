import my_capi as mdl
import numpy as np
import math
import sys

#print(mdl.calc_root(2))
#my_list = [ [1,2,3] , [4,5,6] , [7,8,9]]
#mdl.print_list(my_list)
def print_vector(vector):
    vector_string = ""
    for coordinate in vector:
        vector_string += ("%.4f" % coordinate) + ","
    print(vector_string[:-1])


def print_diagonal_matrix(diag):
    n = len(diag)
    for i in range(n):
        to_print = ('0,' * i + ("%.4f," % diag[i]) + '0,' * (n-i-1))[:-1]
        print(to_print)



def print_matrix(matrix):
    for i in range(len(matrix)):
        print_vector(matrix[i])


def get_vector(file):
    vector = file.readline().split(',')
    if not vector or vector == ['']:
        return -1
    return [float(coordinate) for coordinate in vector]


def initialize_vectors(input_file):
    vectors = []
    with open(input_file) as f:
        while True:
            vector = get_vector(f)
            if vector == -1:
                break
            vectors.append(vector)
    return vectors
    
def average_matrix(mat):
    #gets a matrix (list of lists) and returns the average of all entries in it
    n = len(mat)
    m = len(mat[0])
    total = sum(sum(row) for row in mat)
    return total / (n * m)

def initialize_H(mat_avg, n ,k):
    ans = []
    for i in range(n):
        row = []
        for j in range(k):
            row.append(np.random.uniform(low = 0, high = 2*math.sqrt(mat_avg/k)))
        ans.append(row)
    return ans

def main(arguments):
    if (len(arguments) != 4):
        print("Usage: symnmf.py <K> <goal> <file_name>")
        return -1
    np.random.seed(0)
    k = int(arguments[1])
    goal = arguments[2]
    file_name = arguments[3]
    vectors = initialize_vectors(file_name)
    if goal == 'symnmf': 
        w_mat = mdl.norm(vectors)
        mat_avg = average_matrix(w_mat)
        n = len(vectors)
        h_mat = initialize_H(mat_avg, n ,k)
        print_matrix(mdl.symnmf(w_mat, h_mat, k))
    elif goal == 'sym':
        print_matrix(mdl.sym(vectors))
    elif goal == 'ddg':
        print_diagonal_matrix(mdl.ddg(vectors))
    elif goal == 'norm':
        print_matrix(mdl.norm(vectors))
    else:
        print("An Error Has Occured")
        return -1

if __name__ == "__main__":
    main(sys.argv)

