#ifndef DEFS_H
#define DEFS_H

typedef struct node_v{
    double* vec;
    int vec_length;
    struct node_v* next;
} node_vec;

typedef struct node_d{
    double val;
    struct node_d* next;
} node_double;

typedef struct node_m{
    node_double* seq;
    struct node_m* next;
} node_mat;

#define error_msg "An error has occured\n"

#define MAX_ITER (300)

#define EPSILON (1e-4)

#define BETA (0.5)

#define FREE_LIST(head) free_matrix(head)

void free_matrix(node_mat* head);

#endif 
