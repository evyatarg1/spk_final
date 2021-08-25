#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include "spkmeans.c"
#include <math.h>
#include <Python.h>

static PyObject* c_code(PyObject *self, PyObject *args){
    PyObject *filename, *str_filename, *goal, *str_goal, *res;
    int k;
    char* final_filename, *final_goal;
    Matrix* matrix;

    if (!PyArg_ParseTuple(args, "iOO", &k, &goal, &filename)){
        return NULL;
    }
    str_filename = PyObject_Repr(filename);
    final_filename = PyString_AsString(str_filename);

    str_goal = PyObject_Repr(goal);
    final_goal = PyString_AsString(str_goal);

    matrix = main_logic(k, final_goal, final_filename, 1);

    res = convert_2d_arr_to_python_2d_list(matrix->vertices, matrix->cols, matrix->cols);
    return res;
}

static void


static PyObject* fit(PyObject *self, PyObject *args){
    PyObject *p_centroids, *p_observations, *item, *lst, *res, *current_lst;
    int k, n, dim, max_iter, i, j, k_item;
    double d_item;
    double **observations;
    double **final_centroids;
    int *ini_centroid_indices;

    if (!PyArg_ParseTuple(args, "OOiiii", &p_centroids, &p_observations, &k, &n, &dim, &max_iter)){
        return NULL;
    }
    if (!PyList_Check(p_centroids)) return NULL;
    if (!PyList_Check(p_observations)) return NULL;

    /*observations = a C type 2d array containing all observations*/
    observations = (double**) calloc(n, sizeof(double*));
    assert(observations!=NULL);
    for (i=0; i<n; i++){
        lst = PyList_GetItem(p_observations, i);
        observations[i] = (double*) calloc(dim, sizeof(double));
        assert(observations[i] != NULL);
        for (j = 0; j<dim; j++){
            item = PyList_GetItem(lst, j);
            d_item = PyFloat_AsDouble(item);
            observations[i][j] = d_item;
        }
    }
    /*ini_centroid_indices = a C type 1d array containing indices of initial centroids*/
    ini_centroid_indices = calloc(k, sizeof(int));
    if (ini_centroid_indices == NULL){
        free(observations);
        printf("Cannot allocate memory for ini_centroid_indices\n");
        exit(1);
    }
    for(i=0;i<k;i++){
        item = PyList_GetItem(p_centroids, i);
        k_item = PyLong_AsLong(item);
        ini_centroid_indices[i] = k_item;
    }

    /*final_centroids - a C type 2d array containing final centroids*/
    final_centroids = kmeans(ini_centroid_indices, observations, k, n, dim, max_iter);
    if (final_centroids == NULL){
        return NULL;
    }
    /*convert final_centroids back to python 2d list*/
    res = PyList_New(n);
    for (i=0;i<k;i++){
        current_lst = PyList_New(dim);
        for (j=0;j<dim;j++){
            item = Py_BuildValue("d", final_centroids[i][j]);
            PyList_SetItem(current_lst, j, item);
        }
        PyList_SetItem(res, i, current_lst);
    }

    free(observations);
    free(ini_centroid_indices);

    free(final_centroids);

    return res;

}

Pyobject* convert_2d_arr_to_python_2d_list(double ** arr, int k, int dim){
    Pyobject *res, *current_lst, *item;
    int i, j;

    res = PyList_New(n);
    for (i=0;i<k;i++){
        current_lst = PyList_New(dim);
        for (j=0;j<dim;j++){
            item = Py_BuildValue("d", arr[i][j]);
            PyList_SetItem(current_lst, j, item);
        }
        PyList_SetItem(res, i, current_lst);
    }
}



#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
        FUNC(METH_VARARGS, fit, "c part of code"),
        {NULL, NULL, 0, NULL}   /* sentinel */
};

static struct PyModuleDef _moduledef = {
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        NULL,
        -1,
        _methods
};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void){
    return PyModule_Create(&_moduledef);
}