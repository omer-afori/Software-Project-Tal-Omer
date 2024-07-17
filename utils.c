#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "utils.h"


void free_matrix(double** head, int length){
    int i;
    for(i = 0; i < length; i++){
        free(head[i]);
    }
    free(head);
}

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

void print_matrix(double** head, int length, int width){
    int i;
    for(i = 0; i < length; i++){
        print_vec(head[i], width);
    }
}

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

double inner_product(double* u, double* v, int length){
    int i;
    double sum = 0;
    for(i = 0; i < length; i++){
        sum += u[i] * v[i];
    }
    return sum;
}

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

void print_vec(double* vec, int length){
    int i;
    for(i = 0; i < length; i++){
        printf("%.4f ", vec[i]);
    }
    printf("\n");
}

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

double calc_iteration_formula(double h, double numerator, double denominator){
    return h * (1 - BETA + (BETA*(numerator/denominator)));
}

void calc_iteration(double** curr, double** prev, double** numerator, double** denominator, int length, int width){
    int i, j;
    for (i = 0; i < length; i++){
        for(j = 0; j < width; j++){
            curr[i][j] = calc_iteration_formula(prev[i][j], numerator[i][j], denominator[i][j]);
        }
    }
}

double distance_sqr(double* u, double* v, int length){
    double sum = 0;
    int i;
    for(i = 0; i < length; i++){
        sum += (u[i] - v[i]) * (u[i] - v[i]);
    }
    return sum;
}