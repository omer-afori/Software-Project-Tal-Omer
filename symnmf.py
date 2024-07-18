import symnmfmodule as mdl
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
    

#gets a matrix (list of lists) and returns the average of all entries in it
def average_matrix(mat):
    n = len(mat)
    m = len(mat[0])
    total = sum(sum(row) for row in mat)
    return total / (n * m)


# a function to generate initial H matrix for symmetric NMF - needs to take into account the average of matrix entries and the desired number of components k
def initialize_H(mat_avg, n ,k):
    ans = []
    for i in range(n):
        row = []
        for j in range(k):
            row.append(np.random.uniform(low = 0, high = 2*math.sqrt(mat_avg/k)))
        ans.append(row)
    return ans


# a wrapper function for goal=symnmf - needs to get normalized similarity, and generate initial H
def calculate_symnmf(vectors, k, print_flag):
    w_mat = mdl.norm(vectors)
    mat_avg = average_matrix(w_mat)
    n = len(vectors)
    h_mat = initialize_H(mat_avg, n ,k)
    answer = mdl.symnmf(w_mat, h_mat, k)
    if print_flag == True:
        print_matrix(answer)
    else:
        return answer



def main(arguments):
    if (len(arguments) != 4):
        print("An Error Has Occured")
        return -1
    np.random.seed(0)
    k = int(arguments[1])
    goal = arguments[2]
    file_name = arguments[3]
    vectors = initialize_vectors(file_name)
    if goal == 'symnmf': 
        calculate_symnmf(vectors, k, True)
    elif goal == 'sym':
        print_matrix(mdl.sym(vectors))
    elif goal == 'ddg':
        print_matrix(mdl.ddg(vectors))
    elif goal == 'norm':
        print_matrix(mdl.norm(vectors))
    else:
        print("An Error Has Occured")
        return -1

if __name__ == "__main__":
    main(sys.argv)

