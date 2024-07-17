#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "utils.h"

/**
 * frees an allocated buffer of buffers (matrix)
 * 
 * head: the pointer to the matrix
 * length: the number of rows in the matrix
 */
void free_matrix(double** head, int length){
    int i;
    for(i = 0; i < length; i++){
        free(head[i]);
    }
    free(head);
}

/**
 * allocate and return a matrix, initialized to zero in all coordinates
 * 
 * length: number of pointers in the matrix
 * width: size of each of the buffers in the matrix
 * 
 * returns: the initialized matrix as a double** or NULL for error
 */
double** create_matrix(int length, int width){
    int i, j;
    double** matrix = (double**)malloc(length * sizeof(double*));
    if(matrix == NULL){
        return NULL;
    }
    for(i = 0; i < length; i++){
        matrix[i] = (double*)malloc(width * sizeof(double));
        if(matrix[i] == NULL){
            free_matrix(matrix, i);
            return NULL;
        }
        for(j = 0; j < width; j++){
            matrix[i][j] = 0;
        }
    }
    return matrix;
}

/**
 * prints a vector of doubles, separated by commas. for example:
 * 1.1,2.2,3.3
 * 
 * vec: the vector to print
 * length: its dimension
 */
void print_vec(double* vec, int length){
    int i;
    for(i = 0; i < length-1; i++){
        printf("%.4f,", vec[i]);
    }
    printf("%.4f", vec[length-1]);
    printf("\n");
}

/**
 * prints a matrix in the expected format. for example:
 * 1.1,2.2,3.3
 * 4.4,5.5,6.6
 * 7.7,8.8,9.9
 * 
 * head: the matrix to print
 * length: the number of rows in the matrix
 * width: the number of columns in the matrix
 */
void print_matrix(double** head, int length, int width){
    int i;
    for(i = 0; i < length; i++){
        print_vec(head[i], width);
    }
}

/**
 * implements matrix multiplication (assuming they are of legal dimensions)
 * 
 * k1: length of m1
 * k2: width of m1 and length of m2
 * k3: width of m2
 * 
 * returns: the result of the multiplication of m1 and m2 (a newly allocated matrix!) or NULL for error
 */
double** mat_mult(double** m1, double** m2, int k1, int k2, int k3){
    int i, j, l;
    double sum;
    double** result = create_matrix(k1, k3);
    if(result == NULL){
        return NULL;
    }
    for(i = 0; i < k1; i++){
        for(j = 0; j < k3; j++){
            sum = 0;
            for(l = 0; l < k2; l++){
                sum += m1[i][l] * m2[l][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

/**
 * standard euclidean inner product of vectors of the same dimension
 * 
 * u: first vector for product
 * v: second vector for product
 * length: dimension of the vectors
 * 
 * returns: inner product of u and v
 */
double inner_product(double* u, double* v, int length){
    int i;
    double sum = 0;
    for(i = 0; i < length; i++){
        sum += u[i] * v[i];
    }
    return sum;
}

/**
 * calculates the product of a matrix by its transpose 
 * 
 * mat: the matrix for the operation
 * length: number of rows in mat
 * width: number of columns in mat
 * 
 * returns: a newly allocated matrix containing mat*(mat^T) using inner product, a length x length matrix, or NULL for error
 */
double** mat_mult_transpose(double** mat, int length, int width){
    int i, j;
    double** result = create_matrix(length, length);
    if (result == NULL){
        return NULL;
    }
    for(i = 0; i < length; i++){
        for(j = 0; j < length; j++){
            result[i][j] = inner_product(mat[i], mat[j], width);
        }
    }
    return result;
}

/**
 * counts the number of appearances of a certain char in a string
 * 
 * buf: the string to search in 
 * len: the length of the string
 * key: the char to look for
 * 
 * returns: the number of appearance of the key in buf
 */
int count_char(char* buf, int len, char key){
    int count = 0; 
    int i;
    for(i = 0; i < len; i++){
        if (buf[i] == key){
            count++;
        }
    }
    return count;
}

/**
 * returns the distance squared between two matrices using the standard norm shown in the instructions
 * assumes they are of equal dimensions
 * 
 * matrix_1: first matrix for distance
 * matrix_2: second matrix for distance
 * length: number of rows in the matrices
 * width: number of columns in the matrices
 * 
 * returns: distance squared between matrix_1 and matrix_2
 */
double matrix_distance(double** matrix_1, double** matrix_2, int length, int width){
    double distance = 0;
    double diff;
    int i, j;
    for(i = 0; i < length; i++){
        for(j = 0; j < width; j++){
            diff = matrix_1[i][j] - matrix_2[i][j];
            distance += diff * diff;
        }
    }
    return distance;
}

/**
 * calcualtes a pointwise operation of the algorithm for SYMNMF
 */
double calc_iteration_formula(double h, double numerator, double denominator){
    return h * (1 - BETA + (BETA*(numerator/denominator)));
}

/**
 * calculates a single iteration of the algorithm for factorization as shown in the instructions
 * 
 * curr: the new iteration matrix, assumed to be initialized to zero in all coordinates
 * prev: the previous iteration matrix
 * numerator: the matrix to be in the numerator of the formula
 * denominator: the matrix to be in the denominator of the formula
 * length: the number of rows in the matrices
 * width: the number of columns in the matrices
 */
void calc_iteration(double** curr, double** prev, double** numerator, double** denominator, int length, int width){
    int i, j;
    for (i = 0; i < length; i++){
        for(j = 0; j < width; j++){
            curr[i][j] = calc_iteration_formula(prev[i][j], numerator[i][j], denominator[i][j]);
        }
    }
}

/**
 * calculates the distance squared between two vectors
 * 
 * u: the first vector
 * v: the second vector
 * length: the dimension of the vectors
 * 
 * returns: the distance squared between u and v
 */
double distance_sqr(double* u, double* v, int length){
    double sum = 0;
    int i;
    for(i = 0; i < length; i++){
        sum += (u[i] - v[i]) * (u[i] - v[i]);
    }
    return sum;
}