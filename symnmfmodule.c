#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include "defs.h"
#include "symnmf.h"
#include "utils.h"

/**
 * builds a PyObject from a matrix
 * 
 * mat: the matrix to build from
 * rows: the number of rows in the matrix
 * columns: the number of columns in the matrix
 * 
 * returns: a PyObject pointer for the matrix
 */
static PyObject* python_matrix_from_c(double** mat, int rows, int columns){
    PyObject* result;
    PyObject* row;
    PyObject* item;
    double x;
    int j;
    int k;
    result = PyList_New(rows);
    for (j = 0; j < rows; j++){
        row = PyList_New(columns);
        for (k = 0; k < columns; k++){
            x = mat[j][k];
            item = PyFloat_FromDouble(x);
            PyList_SET_ITEM(row, k, item);
        }
        PyList_SET_ITEM(result, j, row);
    }
    return result;
}

/**
 * builds a matrix from a PyObject containing a list of lists
 * 
 * matrix: the list of lists
 * length: the number of lists
 * width: the number of elements in each list
 * 
 * returns: the matrix as a double**
 */
static double** parse_input_matrix(PyObject *matrix, int length, int width){
    PyObject *item;
    PyObject *inner_item;
    double num;
    int i, j;
    double** result;
    result = create_matrix(length, width);
    if (result == NULL){
        printf(error_msg);
        exit(-1);
    }
    for (i= 0; i < length; i++){
        item = PyList_GetItem(matrix, i);
        for (j= 0; j < width; j++){
            inner_item = PyList_GetItem(item, j);
            num = PyFloat_AsDouble(inner_item);
            result[i][j] = num;
        }
    }
    return result;
}

/**
 * the module function for python: sym()
 * 
 * assumes that the arguments passed from python is a single list of lists, containing the vectors in the database
 * 
 * returns: a PyObject representing the similarity matrix
 */

static PyObject* sym(PyObject* self, PyObject *args)
{
    PyObject *lst;
    PyObject* return_matrix;
    int n;
    int d;
    double **head;
    double **mat;
    if (!PyArg_ParseTuple(args, "O", &lst)){
        printf(error_msg);
        exit(-1);
    }
    n = PyObject_Length(lst);
    if (n <= 0){
        printf(error_msg);
        exit(-1);
    }
    /*d = PyList_Size(PyList_GetItem(lst, 0));*/
    d = PyObject_Length(PyList_GetItem(lst, 0));
    head = parse_input_matrix(lst, n, d);
    mat = sym_wrapper(head, 1, n, d);
    if (mat == NULL){
        free_matrix(head, n);
        printf(error_msg);
        exit(-1);
    }
    return_matrix = python_matrix_from_c(mat, n, n);
    free_matrix(head, n);
    free_matrix(mat, n);
    return return_matrix;
}

/**
 * the module function for python: ddg()
 * 
 * assumes that the arguments passed from python is a single list of lists, containing the vectors in the database
 * 
 * returns: a PyObject representing the diagonal degree matrix
 */
static PyObject* ddg(PyObject* self, PyObject *args)
{
    PyObject *lst;
    PyObject* return_matrix;
    double** head;
    double** diag;
    int n;
    int d;
    if (!PyArg_ParseTuple(args, "O", &lst)){
        printf(error_msg);
        exit(-1);
    }

    n = PyObject_Length(lst);
    if (n <= 0){
        printf(error_msg);
        exit(-1);
    }
    d = PyObject_Length(PyList_GetItem(lst, 0));
    head = parse_input_matrix(lst, n, d);
    diag = ddg_wrapper(head, 1, n, d);
    if (diag == NULL){
        free_matrix(head, n);
        printf(error_msg);
        exit(-1);
    }
    return_matrix = python_matrix_from_c(diag, n, n);
    free_matrix(head, n);
    free_matrix(diag, n);
    return return_matrix;
}

/**
 * the module function for python: norm()
 * 
 * assumes that the arguments passed from python is a single list of lists, containing the vectors in the database
 * 
 * returns: a PyObject representing the normalized similarity matrix
 */
static PyObject* norm(PyObject* self, PyObject *args)
{
    PyObject *lst;
    PyObject *return_mat;
    double** head;
    double** norm_mat;
    int n;
    int d;
    if (!PyArg_ParseTuple(args, "O", &lst)){
        printf(error_msg);
        exit(-1);
    }
    n = PyObject_Length(lst);
    if (n <= 0){
        printf(error_msg);
        exit(-1);
    }
    d = PyObject_Length(PyList_GetItem(lst, 0));

    head = parse_input_matrix(lst, n, d);
    norm_mat = norm_wrapper(head, 1, n, d);
    if (norm_mat == NULL){
        free_matrix(head, n);
        printf(error_msg);
        exit(-1);
    }
    return_mat = python_matrix_from_c(norm_mat, n ,n);
    free_matrix(norm_mat, n);
    free_matrix(head, n);
    return return_mat;
}

/**
 * the module function for python: symnmf()
 * 
 * assumes that the arguments passed from python are the normalized similarity matrix, an initial H, and k 
 * k represents the number of columns in the initial H
 * 
 * returns: a PyObject representing the result of the SYMNMF algorithm; the factorization H
 */
static PyObject* symnmf(PyObject* self, PyObject *args){
    PyObject *W_matrix;
    PyObject *initial_H;
    PyObject *return_mat;
    double** W_matrix_for_C;
    double** initial_H_for_C;
    double** answer_mat;
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
        free_matrix(W_matrix_for_C, n);
        printf(error_msg);
        exit(-1);
    }
    answer_mat = calc_symnmf(W_matrix_for_C, initial_H_for_C, k);
    if (answer_mat == NULL){
        free_matrix(W_matrix_for_C, n);
        free_matrix(initial_H_for_C, n);
        printf(error_msg);
        exit(-1);
    }
    return_mat = python_matrix_from_c(answer_mat, n , k);
    free_matrix(answer_mat, n);
    free_matrix(W_matrix_for_C, n);
    free_matrix(initial_H_for_C, n);
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
    "symnmfmodule",
    NULL,
    -1,
    myMethods
};

PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    PyObject *m;
    m = PyModule_Create(&my_module);
    if (!m) {
        exit(-1);
    }
    return m;
}

