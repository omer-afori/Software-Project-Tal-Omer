#ifndef UTILS_H
#define UTILS_H

#include "defs.h"


void print_vec(double* vec, int length);

void print_matrix(node_mat* matrix);

void print_list_vec(node_vec* head);

void print_list_double(node_double* head);

void print_matrix(node_mat* head);

void print_diagonal(node_double* head, int length);

node_mat* create_matrix(int length, int width);

node_mat* mat_mult(node_mat* m1, node_mat* m2, int k1, int k2, int k3);

void free_list_vec(node_vec* head);

void free_list_double(node_double* head);

void free_matrix(node_mat* head);

int count_char(char* buf, int len, char key);

double get_location(node_mat* matrix, int i, int j);

double matrix_distance(node_mat* matrix_1, node_mat* matrix_2, int length, int width);

node_mat* mat_mult_transpose(node_mat* mat, int length, int width);

double inner_product(node_double* u, node_double* v, int length);

void calc_iteration(node_mat* curr, node_mat* prev, node_mat* numerator, node_mat* denominator, int length, int width);

double calc_iteration_formula(double h, double numerator, double denominator);

#endif 
