#define  _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "symnmf.h"
#include "defs.h"
#include "utils.h"


/**
 * n is the number of vectors in the file
 * d is the dimension of the vectors, which is assumed to be standard across all the vectors
 */
int n = 0;
int d = 0;

/**
 * gets an open file to read from, reads a single line from it, and returns a double* that holds the vector represented in the line 
 * 
 * file: the file to read from
 * 
 * returns: the vector represented in the line as a double*
 */
double* get_vector(FILE *file){
	char *lineptr = NULL;
	size_t size = 0;
	int linelen;
	int my_index = 0;
	double* vector;
	int counter = 0;
    int i;
	linelen = getline(&lineptr, &size, file);
	if(linelen < 0){
		free(lineptr);
		printf(error_msg);
		return NULL;
	}
    if (strncmp(lineptr, "\n", 1) == 0){
        return (double*)1;
    }
	for (i = 0; i < linelen; i++){
		if (lineptr[i] == ',' || lineptr[i] == '\n'){
			lineptr[i] = '\0';
		}
	}

	vector = malloc(sizeof(double)*d);
	if (!vector){
        free(lineptr);
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
    free(lineptr);
	return vector;
}

/**
 * handles general initialization for parameters in the algorithm - the list of vectors (matrix) and the globals n and d
 * 
 * filename: the name of file to read from
 * 
 * returns: double** a buffer of vectors (matrix, n by d)
 */
double** initialize(char* filename){
    double **head;
    initialize_globals_file(filename);
    head = initialize_vectors(filename);
    return head;
}

/**
 * initializes the buffer of vectors from the file
 * 
 * filename: the file to read from
 * 
 * returns: a double** containing all the vectors in the file
 */
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
        if (vector == (double*)1){
            printf("TEST1\n");
            n = i;
            break;
        }
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

/**
 * initializes the globals n and d as explained above
 * 
 * filename: the file to read from
 */
void initialize_globals_file(char *filename){
    size_t size = 0; 
    int count;
    char* line = NULL; 
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

/**
 * calculates a coordinate in the symetric matrix, as instructed in the algorithm
 * 
 * u: the first vector in the formula
 * v: the second vector in the formula
 */
double sym_mat_formula(double* u, double* v){
    double inner;
    if (u == v)
        return 0;
    inner = distance_sqr(u, v, d) / -2.0;
    return exp(inner);
}

/**
 * calculates the first step of the algorithm - similarity matrix
 * 
 * vectors: the list of vectors from the database
 * 
 * returns: the similarity matrix 
 */
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

/**
 * calculates the second step of the algorithm - the diagonal degree matrix
 * 
 * sym_matrix - the similarity matrix from the first step
 * 
 * returns: the diagonal degree matrix
 */
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

/**
 * calculates the (-1/2) power of a diagonal matrix (inplace)
 * 
 * head: the pointer to the matrix
 */
void inverse_sqrt_diag_mat(double** head){
    int i; 
    for(i = 0; i < n; i++){
        head[i][i] = 1.0 / sqrt(head[i][i]);
    }

}

/**
 * calculate the third step of the algorithm - the normalized similarity matrix
 * 
 * sym_matrix: the similarity matrix from the first step of the algorithm
 * ddg_matrix: the diagonal degree matrix from the second step of the algorithm
 * 
 * returns: the normalized similarity matrix W
 */
double** calc_norm(double** sym_matrix, double** ddg_matrix){
    double **DA;
    double **result;
    inverse_sqrt_diag_mat(ddg_matrix);
    DA = mat_mult(ddg_matrix, sym_matrix, n, n, n);
    result = mat_mult(DA, ddg_matrix, n, n, n);
    free_matrix(DA, n);
    return result;
}

/**
 * calculate the final step of the algorithm - the optimized factorization matrix H
 * 
 * w_matrix: the normalized similarity matrix
 * H_matrix: the initial H matrix
 * k: the number of columns in H
 * 
 * returns: the factorization of w_matrix, meaning the optimal H by the method in the instructions
 */
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

/**
 * wrapper function for goal=norm
 * handles the flow of the algorithm - calls the previous steps, returns or prints the normalized similarity matrix
 * 
 * head: the list of vectors in database
 * return_mat: a flag indicating whether to return the matrix, or print it and free it 
 * length: the number of vectors in head
 * dimension: the dimension of each vector in head
 * 
 * returns: 
 *  if return_mat is 1, returns the normalized similar matrix
 *  for an error, returns NULL
 *  returns 1 if return_mat is 0 and there was no error
 */
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

/**
 * wrapper function for goal=ddg
 * handles the flow of the algorithm - calls the previous step and returns or prints the diagonal degree matrix
 * 
 * head: the list of vectors in database
 * return_mat: a flag indicating whether to return the matrix, or print it and free it 
 * length: the number of vectors in head
 * dimension: the dimension of each vector in head
 * 
 * returns: 
 *  if return_mat is 1, returns the diagonal degree matrix
 *  for an error, returns NULL
 *  returns 1 if return_mat is 0 and there was no error
 */
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

/**
 * wrapper function for goal=sym
 * handles the flow of the algorithm and returns or prints the similarity matrix
 * 
 * head: the list of vectors in database
 * return_mat: a flag indicating whether to return the matrix, or print it and free it 
 * length: the number of vectors in head
 * dimension: the dimension of each vector in head
 * 
 * returns: 
 *  if return_mat is 1, returns the similarity matrix
 *  for an error, returns NULL
 *  returns 1 if return_mat is 0 and there was no error
 */
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
    /**
     * in all following cases - the wrapperes return NULL for an error and thus an error message must be printed
     */
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



