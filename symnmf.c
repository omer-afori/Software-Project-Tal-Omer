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
		free(lineptr);
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
        free(lineptr);
		printf(error_msg);
		return NULL;
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

void initialize_globals_node(node_vec* head){
    node_vec* current = head;
    int counter = 0;
    while(current != NULL){
        counter++;
        current = current->next;
    }
    n = counter;
    d = head->vec_length;
    return;
}


node_vec* initialize(char* filename){
    node_vec *head;
    initialize_globals_file(filename);
    head = initialize_vectors(filename);
    return head;
}

node_vec* initialize_vectors(char* filename){
    node_vec* head;
	double* vector;
	node_vec* current_node;
    int i;
    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL){
        printf(error_msg);
        exit(-1);
    }
	head = (node_vec*)malloc(sizeof(node_vec));
    if (!head){
        printf(error_msg);
        exit(-1);
    }
	current_node = head;
    for (i = 0; i  < n; i++){
        if((vector = get_vector(file)) == NULL){
            free_list_vec(head);
            exit(-1);
        }
		current_node->vec = vector;
        current_node->vec_length = d;
        if (i < n - 1){
            current_node->next = (node_vec*)malloc(sizeof(node_vec));
            if (!current_node->next){
                printf(error_msg);
                free_list_vec(head);
                exit(-1);
            }
        }
        else {
            current_node->next = NULL;
        }
		current_node = current_node->next;
    }
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

double distance_sqr(double* u, double* v){
    double sum = 0;
    int i;
    for(i = 0; i < d; i++){
        sum += (u[i] - v[i]) * (u[i] - v[i]);
    }
    return sum;
}


double sym_mat_formula(double* u, double* v){
    double inner;
    if (u == v)
        return 0;
    inner = distance_sqr(u, v) / -2.0;
    return exp(inner);
}

node_mat* calc_sym(node_vec* vectors){
    float calc;
    int i, j; 
    node_double* row_iter;
    node_vec* vectors_i = vectors;
    node_vec* vectors_j = vectors;
    node_mat* head;
    node_mat* curr;
    head = create_matrix(n, n);
    curr = head;
    if (head == NULL){
        return NULL;
    }
    for(i = 0; i < n; i++){
        row_iter = curr->seq;
        vectors_j = vectors;
        for(j = 0; j < n; j++){
            calc = sym_mat_formula(vectors_i->vec, vectors_j->vec);
            row_iter->val = calc;
            row_iter = row_iter->next;
            vectors_j = vectors_j->next;
        }
        vectors_i = vectors_i->next;
        curr = curr->next;
    }
    return head;
}

node_double* calc_ddg(node_mat* sym_matrix){
    node_double* head;
    node_double* tmp;
    node_mat* column_iter = sym_matrix;
    node_double* row_iter;
    int i, j;
    double sum;
    head = malloc(sizeof(node_double));
    tmp = head;
    if (head == NULL){
        printf(error_msg);
        return NULL;
    }
    for (i = 0; i < n; i++){
        sum = 0;
        row_iter = column_iter->seq;
        for (j = 0; j < n; j++){
            sum += row_iter->val;
            row_iter = row_iter->next;
        }
        column_iter = column_iter->next;
        tmp->val = sum;
        if (i < n - 1){
            tmp->next = malloc(sizeof(node_double));
            if (tmp->next == NULL){
                printf(error_msg);
                free_list_double(head);
                return NULL;
            }
        }
        else {
            tmp->next = NULL;
        }
        tmp = tmp->next;
    }
    return head;
}

void inverse_sqrt_diag_mat(node_double* head){
    int i; 
    for(i = 0; i < n; i++){
        head->val = 1.0 / sqrt(head->val);
        head = head->next;
    }
}

node_mat* calc_norm(node_mat* sym_matrix, node_double* ddg_matrix){
    int i, j;
    node_double* d_i = ddg_matrix;
    node_double* d_j;
    node_mat* iter = sym_matrix;
    node_double* row_iter;
    inverse_sqrt_diag_mat(ddg_matrix);
    for(i = 0; i < n; i++){
        row_iter = iter->seq;
        d_j = ddg_matrix;
        for(j = 0; j < n; j++){
            row_iter->val = row_iter->val * d_i->val * d_j->val;
            row_iter = row_iter->next;
            d_j = d_j->next;
        }
        d_i = d_i->next;
        iter = iter->next;
    }
    return sym_matrix;
}

node_mat* calc_symnmf(node_mat* w_matrix, node_mat* H_matrix, int k){
    int iter = 0; 
    node_mat* h_t = H_matrix;
    node_mat *h_t_1;
    node_mat *neumerator, *denominator, *tmp;
    double dist;
    do{
        neumerator = mat_mult(w_matrix, h_t, n, n, k);
        if (neumerator == NULL){
            return NULL;
        }
        tmp = mat_mult_transpose(h_t, n, k);
        if (tmp == NULL){
            free_matrix(neumerator);
            return NULL;
        }
        denominator = mat_mult(tmp, h_t, n, n, k);
        if (denominator == NULL){
            free_matrix(neumerator);
            free_matrix(tmp);
            return NULL;
        }
        h_t_1 = create_matrix(n, k);
        if (h_t_1 == NULL){
            free_matrix(neumerator);
            free_matrix(tmp);
            free_matrix(denominator);
            return NULL;
        }

        calc_iteration(h_t_1, h_t, neumerator, denominator, n, k);
        dist = matrix_distance(h_t, h_t_1, n, k);
        if (iter != 0){
            free_matrix(h_t);
        }
        free_matrix(neumerator);
        free_matrix(denominator);
        free_matrix(tmp);
        h_t = h_t_1;
        iter++;
    }while (iter < MAX_ITER && dist >= EPSILON);
    return h_t_1;
}

node_mat* norm_wrapper(node_vec* head, int globals_flag, int return_mat){
    node_mat* sym_matrix;
    node_double* ddg_matrix;
    node_mat* norm_matrix;
    if (globals_flag == 1){
        initialize_globals_node(head);
    }
    sym_matrix = calc_sym(head);
    if (sym_matrix == NULL){
        return NULL;
    }
    ddg_matrix = calc_ddg(sym_matrix);
    if (ddg_matrix == NULL){
        free_matrix(sym_matrix);
        return NULL;
    }
    norm_matrix = calc_norm(sym_matrix, ddg_matrix);
    free_list_double(ddg_matrix);
    if (norm_matrix == NULL){
        return NULL;
    }
    if (return_mat == 1){
        return norm_matrix;
    }
    else{
        print_matrix(norm_matrix);
        free_matrix(norm_matrix);
        return (node_mat*)1;
    }
}

node_double* ddg_wrapper(node_vec* head, int globals_flag, int return_mat){
    node_mat* sym_matrix;
    node_double* ddg_matrix;
    if (globals_flag == 1){
        initialize_globals_node(head);
    }
    sym_matrix = calc_sym(head);
    if (sym_matrix == NULL){
        return NULL;
    }
    ddg_matrix = calc_ddg(sym_matrix);
    free_matrix(sym_matrix);
    if (ddg_matrix == NULL){
        return NULL;
    }
    if (return_mat == 1){
        return ddg_matrix;
    }
    else{
        print_diagonal(ddg_matrix, n);
        free_list_double(ddg_matrix);
        return (node_double*)1;
    }
}

node_mat* sym_wrapper(node_vec* head, int globals_flag, int return_mat){
    node_mat* sym_matrix;
    if (globals_flag == 1){
        initialize_globals_node(head);
    }
    sym_matrix = calc_sym(head);
    if (sym_matrix == NULL){
        return NULL;
    }
    if (return_mat == 1){
        return sym_matrix;
    }
    else{
        print_matrix(sym_matrix);
        free_matrix(sym_matrix);
        return (node_mat*)1;
    }

}

int main(int argc, char *argv[]) {
    char* goal;
    char* filename;
    node_vec* head;

    if (argc!= 3) { 
        printf(error_msg);
        exit(-1);
    }

    goal = argv[1];
    filename = argv[2];
    if (strncmp(goal, "sym", 4) == 0){
        head = initialize(filename);
        if (sym_wrapper(head, 0, 0) == NULL){
            free_list_vec(head);
            printf(error_msg);
            exit(-1);
        }
    }
    else if (strncmp(goal,"ddg", 4) == 0){
        head = initialize(filename);
        if (ddg_wrapper(head, 0, 0) == NULL){
            free_list_vec(head);
            printf(error_msg);
            exit(-1);
        }
    }
    else if (strncmp(goal, "norm", 5) == 0){
        head = initialize(filename);
        if (norm_wrapper(head, 0, 0) == NULL){
            free_list_vec(head);
            printf(error_msg);
            exit(-1);
        }
    }
    else{
        printf(error_msg);
        exit(-1);
    }
    free_list_vec(head);
    return 1;
}



