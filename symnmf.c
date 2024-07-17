#define  _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "symnmf.h"
#include "defs.h"
#include "utils.h"

int n = 0;
int d = 0;

double* get_vector(FILE *file){
	char *lineptr = 0;
	size_t size;
	int linelen;
	int my_index = 0;
	double* vector;
	int counter = 0;
    int i;
	linelen = getline(&lineptr, &size, file);
	if(linelen < 0){
		/*free(lineptr);*/
		printf(error_msg);
		return NULL;
	}

	for (i = 0; i < linelen; i++){
		if (lineptr[i] == ',' || lineptr[i] == '\n'){
			lineptr[i] = '\0';
		}
	}

	vector = malloc(sizeof(double)*d);
	if (!vector){
        /*free(lineptr);*/
		printf(error_msg);
		return NULL;
	}
    for (i = 0; i < d; i++){
        vector[i] = 0.0;
    }
	for(i = 0; i < linelen; i++){
		if (lineptr[i] == '\0'){
			vector[counter] = atof(lineptr + my_index);
			my_index = i+1;
			counter++;
		}
	}
    /*free(lineptr);*/
	return vector;
}


double** initialize(char* filename){
    double **head;
    initialize_globals_file(filename);
    head = initialize_vectors(filename);
    return head;
}

double** initialize_vectors(char* filename){
    double** head;
	double* vector;
    int i;
    int j;
    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL){
        printf(error_msg);
        exit(-1);
    }
	head = (double**)malloc(sizeof(double*)*n);
    if (head == NULL){
        printf(error_msg);
        fclose(file);
        exit(-1);
    }
    for (i=0; i<n; i++){
        vector = get_vector(file);
        if (vector == NULL){
            for (j = 0; j < i; j ++){
                free(head[j]);
            }
            free(head);
            printf(error_msg);
            fclose(file);
            exit(-1);
        }
        head[i] = vector;
    }
    fclose(file);
    return head;
}   

void initialize_globals_file(char *filename){
    size_t size; 
    int count;
    char* line; 
    FILE *file;
    int linelen;
    char ch;
    file = fopen(filename, "r");
    if (file == NULL){
        printf(error_msg);
        exit(-1);
    }
    linelen = getline(&line, &size, file);
    if (linelen < 0){
        printf(error_msg);
        free(line);
        exit(-1);
    }
    d = count_char(line, linelen, ',') + 1;
    free(line);
    count = 1; /* assume exactly one blank line at the end of the file */
    while ((ch = fgetc(file)) != EOF){
        if (ch == '\n'){
            count++;
        }
    }
    n = count;
    fclose(file);
}




double sym_mat_formula(double* u, double* v){
    double inner;
    if (u == v)
        return 0;
    inner = distance_sqr(u, v, d) / -2.0;
    return exp(inner);
}

double** calc_sym(double** vectors){
    int i, j; 
    double** head;
    head = create_matrix(n, n);
    if (head == NULL){
        return NULL;
    }
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            head[i][j] = sym_mat_formula(vectors[i], vectors[j]);
        }
    }
    return head;
}

double** calc_ddg(double** sym_matrix){
    double** head;
    int i, j;
    double sum;
    head = create_matrix(n, n);
    if (head == NULL){
        return NULL;
    }
    for (i = 0; i < n; i++){
        sum = 0.0;
        for (j = 0; j < n; j++){
            sum += sym_matrix[i][j];
        }
        head[i][i] = sum;
    }
    return head;
}

void inverse_sqrt_diag_mat(double** head){
    int i; 
    for(i = 0; i < n; i++){
        head[i][i] = 1.0 / sqrt(head[i][i]);
    }

}

double** calc_norm(double** sym_matrix, double** ddg_matrix){
    double **DA;
    double **result;
    inverse_sqrt_diag_mat(ddg_matrix);
    DA = mat_mult(ddg_matrix, sym_matrix, n, n, n);
    result = mat_mult(DA, ddg_matrix, n, n, n);
    free_matrix(DA, n);
    return result;
}

double** calc_symnmf(double** w_matrix, double** H_matrix, int k){
    int iter = 0; 
    double** h_t = H_matrix;
    double** h_t_1;
    double** neumerator;
    double** denominator;
    double** tmp;
    double dist;
    do{
        neumerator = mat_mult(w_matrix, h_t, n, n, k);
        if (neumerator == NULL){
            return NULL;
        }
        tmp = mat_mult_transpose(h_t, n, k);
        if (tmp == NULL){
            free_matrix(neumerator, n);
            return NULL;
        }
        denominator = mat_mult(tmp, h_t, n, n, k);
        if (denominator == NULL){
            free_matrix(neumerator, n);
            free_matrix(tmp, n);
            return NULL;
        }
        h_t_1 = create_matrix(n, k);
        if (h_t_1 == NULL){
            free_matrix(neumerator, n);
            free_matrix(tmp, n);
            free_matrix(denominator, n);
            return NULL;
        }
        calc_iteration(h_t_1, h_t, neumerator, denominator, n, k);
        dist = matrix_distance(h_t, h_t_1, n, k);
        if (iter != 0){
            free_matrix(h_t, n);
        }
        free_matrix(neumerator, n);
        free_matrix(denominator, n);
        free_matrix(tmp, n);
        h_t = h_t_1;
        iter++;
    }while (iter < MAX_ITER && dist >= EPSILON);
    return h_t_1;
}

double** norm_wrapper(double** head, int return_mat, int length, int dimension){
    double** sym_matrix;
    double** ddg_matrix;
    double** norm_matrix;
    n = length;
    d = dimension;
    sym_matrix = calc_sym(head);
    if (sym_matrix == NULL){
        return NULL;
    }
    ddg_matrix = calc_ddg(sym_matrix);
    if (ddg_matrix == NULL){
        free_matrix(sym_matrix, n);
        return NULL;
    }
    norm_matrix = calc_norm(sym_matrix, ddg_matrix);
    free_matrix(ddg_matrix, n);
    free_matrix(sym_matrix, n);
    if (norm_matrix == NULL){
        return NULL;
    }
    if (return_mat == 1){
        return norm_matrix;
    }
    else{
        print_matrix(norm_matrix, n, n);
        free_matrix(norm_matrix, n);
        return (double**)1;
    }
}

double** ddg_wrapper(double** head, int return_mat, int length, int dimension){
    double** sym_matrix;
    double** ddg_matrix;
    n = length;
    d = dimension;
    sym_matrix = calc_sym(head);
    if (sym_matrix == NULL){
        return NULL;
    }
    ddg_matrix = calc_ddg(sym_matrix);
    free_matrix(sym_matrix, n);
    if (ddg_matrix == NULL){
        return NULL;
    }
    if (return_mat == 1){
        return ddg_matrix;
    }
    else{
        print_matrix(ddg_matrix, n, n);
        free_matrix(ddg_matrix, n);
        return (double**)1;
    }
}

double** sym_wrapper(double** head, int return_mat, int length, int dimension){
    double** sym_matrix;
    n = length;
    d = dimension;
    sym_matrix = calc_sym(head);
    if (sym_matrix == NULL){
        return NULL;
    }
    if (return_mat == 1){
        return sym_matrix;
    }
    else{
        print_matrix(sym_matrix, n, n);
        free_matrix(sym_matrix, n);
        return (double**)1;
    }
}

int main(int argc, char *argv[]) {
    char* goal;
    char* filename;
    double** head;

    if (argc!= 3) { 
        printf(error_msg);
        exit(-1);
    }

    goal = argv[1];
    filename = argv[2];
    if (strncmp(goal, "sym", 4) == 0){
        head = initialize(filename);
        if (sym_wrapper(head, 0, n, d) == NULL){
            free_matrix(head, n);
            printf(error_msg);
            exit(-1);
        }
    }
    else if (strncmp(goal,"ddg", 4) == 0){
        head = initialize(filename);
        if (ddg_wrapper(head, 0, n, d) == NULL){
            free_matrix(head, n);
            printf(error_msg);
            exit(-1);
        }
    }
    else if (strncmp(goal, "norm", 5) == 0){
        head = initialize(filename);
        if (norm_wrapper(head, 0, n, d) == NULL){
            free_matrix(head, n);
            printf(error_msg);
            exit(-1);
        }
    }
    else{
        printf(error_msg);
        exit(-1);
    }
    free_matrix(head, n);
    return 1;
}



