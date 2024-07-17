#ifndef SYM_NMF_H
#define SYM_NMF_H

#include "defs.h"

double distance_dqr(double* u, double* v);

double sym_mat_formula(double* u, double* v);

node_vec* initialize(char* filename);

void initialize_globals_node(node_vec* head);

void initialize_globals_file(char* filename);

node_vec* initialize_vectors(char* filename);

node_mat* calc_sym(node_vec* vectors);

node_double* calc_ddg(node_mat* sym_matrix);

node_mat* calc_norm(node_mat* sym_matrix, node_double* ddg_matrix);

node_mat* calc_symnmf(node_mat* w_mat, node_mat* H_mat, int k);

node_mat* sym_wrapper(node_vec* head, int globals_flag, int return_mat);

node_double* ddg_wrapper(node_vec* head, int globals_flag, int return_mat);

node_mat* norm_wrapper(node_vec* head, int globals_flag, int return_mat);

double* get_vector(FILE *file);

#endif 
