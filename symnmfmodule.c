#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include "defs.h"
#include "symnmf.h"
#include "utils.h"



static PyObject* list_from_node_double(node_double* head, int length){
    PyObject* result;
    PyObject* item;
    double x;
    int j;
    result = PyList_New(length);
    for (j = 0; j < length; j++){
        x = head->val;
        item = PyFloat_FromDouble(x);
        PyList_SET_ITEM(result, j, item);
        head = head->next;
    }
    return result;
}

static PyObject* matrix_from_node_mat(node_mat* mat, int rows, int columns){
    PyObject* result;
    PyObject* row;
    PyObject* item;
    node_double* row_iter;
    double x;
    int j;
    int k;
    result = PyList_New(rows);
    for (j = 0; j < rows; j++){
        row_iter = mat->seq;
        row = PyList_New(columns);
        for (k = 0; k < columns; k++){
            x = row_iter->val;
            item = PyFloat_FromDouble(x);
            PyList_SET_ITEM(row, k, item);
            row_iter = row_iter->next;
        }
        PyList_SET_ITEM(result, j, row);
        mat = mat->next;
    }
    return result;
}

static node_mat* parse_input_matrix(PyObject *matrix, int length, int width){
    node_mat* result;
    node_mat* iterator;
    node_double* row_iter;
    PyObject *item;
    PyObject *inner_item;
    double num;
    int i, j;
    result = create_matrix(length, width);
    iterator = result;
    if (result == NULL){
        printf(error_msg);
        exit(-1);
    }
    for (i= 0; i < length; i++){
        row_iter = iterator->seq;
        item = PyList_GetItem(matrix, i);
        for (j= 0; j < width; j++){
            inner_item = PyList_GetItem(item, j);
            num = PyFloat_AsDouble(inner_item);
            row_iter->val = num;
            row_iter = row_iter->next;
        }
        iterator = iterator->next;
    }
    return result;
}

static node_vec* parse_input_vectors(PyObject *lst, int n){
    node_vec *tmp;
    node_vec *head;
    PyObject *item;
    PyObject *inner_item;
    double num;
    int j;
    int k;
    head = (node_vec*)malloc(sizeof(node_vec));
    if (head == NULL){
        printf(error_msg);
        exit(-1);
    }
    tmp = head;
    for (j = 0; j < n; j++){
        if (j < n - 1){
            tmp->next = (node_vec*)malloc(sizeof(node_vec));
            if (tmp->next == NULL){
                printf(error_msg);
                free_list_vec(head);
                exit(-1);
            }
        }
        else{
            tmp->next = NULL;
        }
        item = PyList_GetItem(lst, j);
        int d = PyObject_Length(item);
        tmp->vec_length = d;
        tmp->vec = (double*)malloc(sizeof(double) * d);
        if (tmp->vec == NULL){
            free(tmp->next);
            tmp->next = NULL;
            free_list_vec(head);
            printf(error_msg);
            exit(-1);
        }
        for (k = 0; k < d; k++) {
            inner_item = PyList_GetItem(item, k);
            num = PyFloat_AsDouble(inner_item);
            tmp->vec[k] = num;
        }
        tmp = tmp->next;
    }
    return head;
}

static PyObject* sym(PyObject* self, PyObject *args)
{
    PyObject *lst;
    PyObject* return_matrix;
    node_vec *head;
    node_mat *mat;
    int n;
    if (!PyArg_ParseTuple(args, "O", &lst)){
        printf(error_msg);
        exit(-1);
    }
    n = PyObject_Length(lst);
    if (n <= 0){
        printf(error_msg);
        exit(-1);
    }
    head = parse_input_vectors(lst, n);
    mat = sym_wrapper(head, 1, 1);
    if (mat == NULL){
        free_list_vec(head);
        printf(error_msg);
        exit(-1);
    }
    return_matrix = matrix_from_node_mat(mat, n, n);
    free_list_vec(head);
    free_matrix(mat);
    return return_matrix;
}

static PyObject* ddg(PyObject* self, PyObject *args)
{
    PyObject *lst;
    node_vec *head;
    PyObject* return_matrix;
    node_double* diag;
    int n;
    if (!PyArg_ParseTuple(args, "O", &lst)){
        printf(error_msg);
        exit(-1);
    }

    n = PyObject_Length(lst);
    if (n <= 0){
        printf(error_msg);
        exit(-1);
    }
    head = parse_input_vectors(lst, n);
    diag = ddg_wrapper(head, 1, 1);
    if (diag == NULL){
        free_list_vec(head);
        printf(error_msg);
        exit(-1);
    }
    return_matrix = list_from_node_double(diag, n);
    free_list_vec(head);
    free_list_double(diag);
    return return_matrix;
}

static PyObject* norm(PyObject* self, PyObject *args)
{
    PyObject *lst;
    PyObject *return_mat;
    node_vec *head;
    node_mat* norm_mat;
    int n;

    if (!PyArg_ParseTuple(args, "O", &lst)){
        printf(error_msg);
        exit(-1);
    }
    n = PyObject_Length(lst);
    if (n <= 0){
        printf(error_msg);
        exit(-1);
    }
    head = parse_input_vectors(lst, n);
    norm_mat = norm_wrapper(head, 1, 1);
    if (norm_mat == NULL){
        free_list_vec(head);
        printf(error_msg);
        exit(-1);
    }
    return_mat = matrix_from_node_mat(norm_mat, n ,n);

    free_matrix(norm_mat);
    free_list_vec(head);
    return return_mat;
}

static PyObject* symnmf(PyObject* self, PyObject *args){
    PyObject *W_matrix;
    PyObject *initial_H;
    PyObject *return_mat;
    node_mat* W_matrix_for_C;
    node_mat* initial_H_for_C;
    node_mat* answer_mat;
    int k;
    int n;
    if (!PyArg_ParseTuple(args, "OOi", &W_matrix, &initial_H, &k)){
        printf(error_msg);
        exit(-1);
    }
    n = PyObject_Length(W_matrix);
    if (n <= 0){
        printf(error_msg);
        exit(-1);
    }
    W_matrix_for_C = parse_input_matrix(W_matrix, n, n);

    if (W_matrix_for_C == NULL){
        printf(error_msg);
        exit(-1);
    }
    
    initial_H_for_C = parse_input_matrix(initial_H, n, k);
    
    if (initial_H_for_C == NULL){
        free_matrix(W_matrix_for_C);
        printf(error_msg);
        exit(-1);
    }
    answer_mat = calc_symnmf(W_matrix_for_C, initial_H_for_C, k);
    if (answer_mat == NULL){
        free_matrix(W_matrix_for_C);
        free_matrix(initial_H_for_C);
        printf(error_msg);
        exit(-1);
    }
    return_mat = matrix_from_node_mat(answer_mat, n , k);
    free_matrix(answer_mat);
    free_matrix(W_matrix_for_C);
    free_matrix(initial_H_for_C);
    return return_mat;
}

static PyMethodDef myMethods[] = {

    {"symnmf",
    (PyCFunction)symnmf,
    METH_VARARGS,
    PyDoc_STR("symnmf")},

    {"norm",
    (PyCFunction)norm,
    METH_VARARGS,
    PyDoc_STR("norm")},

    {"ddg",
    (PyCFunction)ddg,
    METH_VARARGS,
    PyDoc_STR("ddg")},

    {"sym",
    (PyCFunction)sym,
    METH_VARARGS,
    PyDoc_STR("sym")},
    
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef my_module = {
    PyModuleDef_HEAD_INIT, 
    "my_capi",
    NULL,
    -1,
    myMethods
};

PyMODINIT_FUNC PyInit_my_capi(void) {
    PyObject *m;
    m = PyModule_Create(&my_module);
    if (!m) {
        exit(-1);
    }
    return m;
}

