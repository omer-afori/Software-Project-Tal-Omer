#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "utils.h"

void print_list_vec(node_vec *head){
    node_vec *current = head;
    int i;
    while(current!= NULL){
        for (i = 0; i < current->vec_length - 1; i++){
            printf("%.4f,", current->vec[i]);
        }
        printf("%.4f\n", current->vec[current->vec_length - 1]);
        current = current->next;
    }
}

void print_list_double(node_double *head){
    node_double *current = head;
    if (current == NULL){
        printf(error_msg);
        return;
    }
    while(current->next != NULL){
        printf("%.4f,", current->val);
        current = current->next;
    }
    printf("%.4f", current->val);
}


void free_list_vec(node_vec* head){
    node_vec *current = head;
    node_vec *next;
    while(current!= NULL){
        next = current->next;
        if (current->vec)
            free(current->vec);
        free(current);
        current = next;
    }
}

void free_list_double(node_double* head){
    node_double *current = head;
    node_double *next;
    while(current!= NULL){
        next = current->next;
        free(current);
        current = next;
    }
}

void free_matrix(node_mat *head){
    node_mat *current = head;
    node_mat *next;
    while(current != NULL){
        next = current->next;
        if (current->seq)
            free_list_double(current->seq);
        free(current);
        current = next;
    }
}

void print_diagonal(node_double* head, int length){
    node_double *current = head;
    int i, j;
    for(i = 0; i < length; i++){
        for(j = 0; j < length; j++){
            if(j == i){
                printf("%.4f", current->val);
                current = current->next;
            } 
            else {
                printf("0.0000");
            }
            if (j != length - 1){
                printf(",");
            }
        }
        printf("\n");
    }
}

node_mat* create_matrix(int length, int width){
    int i, j; 
    node_double* row_iter;
    node_mat* curr;
    node_mat* head = malloc(sizeof(node_mat));
    if (head == NULL){
        printf(error_msg);
        return NULL;
    }
    curr = head;
    for(i = 0; i < length; i++){
        curr->seq = malloc(sizeof(node_double));
        if (curr->seq == NULL){
            free_matrix(head);
            printf(error_msg);
            return NULL;
        }
        curr->seq->next = malloc(sizeof(node_double));
        if (curr->seq->next == NULL){
            free_matrix(head);
            printf(error_msg);
            return NULL;
        }
        row_iter = curr->seq;
        for(j = 0; j < width; j++){
            row_iter->val = 0;
            if (j < width - 1){
                row_iter->next = malloc(sizeof(node_double));
                if (row_iter->next == NULL){
                    free_matrix(head);
                    printf(error_msg);
                    return NULL;
                }
            }
            else {
                row_iter->next = NULL;
            }
            row_iter = row_iter->next;
        }
        if (i < length - 1){
            curr->next = malloc(sizeof(node_double));
            if (curr->next == NULL){
                free_matrix(head);
                printf(error_msg);
                return NULL;
            } 
        }
        else {
            curr->next = NULL;
        }
        curr = curr->next;
    }
    return head;
}

void print_matrix(node_mat *head){
    node_mat *current = head;
    while(current != NULL){
        if (current->seq){
            print_list_double(current->seq);
        }
        printf("\n");
        current = current->next;
    }
}

double get_location(node_mat* matrix, int i, int j){
    int iter_i;
    int iter_j;
    node_double* row_iter;
    node_mat* mat_iter = matrix;
    for (iter_i = 0; iter_i < i; iter_i++){
        mat_iter = mat_iter->next;
    }
    row_iter = mat_iter->seq;
    for (iter_j = 0; iter_j < j; iter_j++){
        row_iter = row_iter->next;
    }
    return row_iter->val;
}

node_mat* mat_mult(node_mat* m1, node_mat* m2, int k1, int k2, int k3){
    int i, j, l;
    double sum;
    node_mat* result;
    node_mat* iter;
    node_double* row_iter;
    result = create_matrix(k1, k3);
    if (result == NULL){
        return NULL;
    }
    iter = result;
    for (i = 0; i < k1; i++){
        row_iter = iter->seq;
        for (j = 0; j < k3; j++){
            sum = 0;
            for (l = 0; l < k2; l++){
                sum += get_location(m1, i, l) * get_location(m2, l, j);
            }
            row_iter->val = sum;
            row_iter = row_iter->next;
        }
        iter = iter->next;
    }
    return result;
}

double inner_product(node_double* u, node_double* v, int length){
    int i;
    double sum = 0;
    for(i = 0; i < length; i++){
        sum += u->val * v->val;
        u = u->next;
        v = v->next;
    }
    return sum;
}

node_mat* mat_mult_transpose(node_mat* mat, int length, int width){
    node_mat* iter_i = mat;
    node_mat* iter_j = mat;
    node_mat* ans = create_matrix(length, length);
    node_mat* head = ans;
    node_double* ans_row_iter;
    int i, j;
    if (ans == NULL){
        return NULL;
    }
    for (i = 0; i < length; i++){
        iter_j = mat;
        ans_row_iter = head->seq;
        for (j = 0; j < length; j++){
            ans_row_iter->val = inner_product(iter_i->seq, iter_j->seq, width);
            iter_j = iter_j->next;
            ans_row_iter = ans_row_iter->next;
        }
        iter_i = iter_i->next;
        head = head->next;
    }
    return ans;
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


double matrix_distance(node_mat* matrix_1, node_mat* matrix_2, int length, int width){
    /*distance is defined as the sum of the squared distances between coordinates*/
    double distance = 0;
    double diff;
    node_double* row_iter_1;
    node_double* row_iter_2;
    int i, j;
    for (i = 0; i < length; i++){
        row_iter_1 = matrix_1->seq;
        row_iter_2 = matrix_2->seq;
        for (j = 0; j < width; j++){
            diff = row_iter_1->val - row_iter_2->val;
            distance += diff * diff;
            row_iter_1 = row_iter_1->next;
            row_iter_2 = row_iter_2->next;
        }
        matrix_1 = matrix_1->next;
        matrix_2 = matrix_2->next;
    }
    return distance;
}

double calc_iteration_formula(double h, double numerator, double denominator){
    return h * (1 - BETA + (BETA*(numerator/denominator)));
}

void calc_iteration(node_mat* curr, node_mat* prev, node_mat* numerator, node_mat* denominator, int length, int width){
    node_double *curr_iter, *prev_iter, *numerator_iter, *denominator_iter;
    int i, j;
    for (i = 0; i < length; i++){
        curr_iter = curr->seq;
        prev_iter = prev->seq;
        numerator_iter = numerator->seq;
        denominator_iter = denominator->seq;
        for (j = 0; j < width; j++){
            curr_iter->val = calc_iteration_formula(prev_iter->val, numerator_iter->val, denominator_iter->val);
            curr_iter = curr_iter->next;
            prev_iter = prev_iter->next;
            numerator_iter = numerator_iter->next;
            denominator_iter = denominator_iter->next;
        }
        curr = curr->next;
        prev = prev->next;
        numerator = numerator->next;
        denominator = denominator->next;
    }
}