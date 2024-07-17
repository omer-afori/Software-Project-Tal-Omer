#ifndef SYM_NMF_H
#define SYM_NMF_H

#include "defs.h"

double sym_mat_formula(double* u, double* v);

double** initialize(char* filename);

void initialize_globals_file(char* filename);

double** initialize_vectors(char* filename);

double** calc_sym(double** vectors);

double** calc_ddg(double** sym_matrix);

double** calc_norm(double** sym_matrix, double** ddg_matrix);

double** calc_symnmf(double** w_mat, double** H_mat, int k);

double** sym_wrapper(double** head, int return_mat, int length, int dimension);

double** ddg_wrapper(double** head, int return_mat, int length, int dimension);

double** norm_wrapper(double** head, int return_mat, int length, int dimension);

double* get_vector(FILE *file);

#endif 
