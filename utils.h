#ifndef UTILS_H
#define UTILS_H

#include "defs.h"


void print_vec(double* vec, int length);

void print_matrix(double** matrix, int length, int width);

double** create_matrix(int length, int width);

double** mat_mult(double** m1, double** m2, int k1, int k2, int k3);

void free_matrix(double** head, int length);

int count_char(char* buf, int len, char key);

double matrix_distance(double** matrix_1, double** matrix_2, int length, int width);

double** mat_mult_transpose(double** mat, int length, int width);

double inner_product(double* u, double* v, int length);

void calc_iteration(double** curr, double** prev, double** numerator, double** denominator, int length, int width);

double calc_iteration_formula(double h, double numerator, double denominator);

double distance_sqr(double* u, double* v, int length);

#endif 
