#define PY_SSIZE_T_CLEAN

#include "spkmeans.h"
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

static PyObject* c_code(PyObject *self, PyObject *args);
static PyObject* c_code_kmeans(PyObject *self, PyObject *args);


static PyObject* c_code(PyObject *self, PyObject *args){
    PyObject *p_data, *item, *lst, *current_lst, *res;
    int k, goal, n, dim, i, j;
    double** data, d_item;
    Matrix* matrix;
    Point** point_arr;

    printf("got to c\n");
    exit(1);
    if (!PyArg_ParseTuple(args, "iOiii", &k, &data, &goal, &n, &dim)){
        return NULL;
    }

    if (!PyList_Check(data)) return NULL;

    /*convert python list to c 2d array */
    data = (double **) calloc(n, sizeof (double *));
    assert(data !=NULL);

    for (i=0; i<n; i++){
        lst = PyList_GetItem(p_data, i);
        data[i] = (double*) calloc(dim, sizeof(double));
        assert(data[i] != NULL);
        for (j = 0; j<dim; j++){
            item = PyList_GetItem(lst, j);
            d_item = PyFloat_AsDouble(item);
            data[i][j] = d_item;
        }
    }

    /*covert c d2 array to point array*/
    point_arr = (Point**) calloc(n, sizeof(Point*));
    assert(point_arr != NULL);
    for (i=0; i<n; i++){
        point_arr[i] = (Point*) calloc(1, sizeof (Point));
        assert(point_arr[i] != NULL);
        point_arr[i]->coordinates = calloc(dim+1, sizeof(double));
        assert(point_arr[i]->coordinates != NULL);
        for (j=0; j<dim; j++){
            point_arr[i]->coordinates[j] = data[i][j];
        }
    }

    if (goal == 0){
        matrix = main_logic(k, "spk", point_arr, n, dim, 1);

        for (i=0;i<k;i++) {
            current_lst = PyList_New(dim);
            for (j = 0; j < dim; j++) {
                item = Py_BuildValue("d", matrix->vertices[i][j]);
                PyList_SetItem(current_lst, j, item);
            }
            PyList_SetItem(res, i, current_lst);
        }

        return res;
    }
    else if (goal == 1){
        main_logic(k, "wam", point_arr, n, dim, 1);
    }
    else if (goal == 2){
        main_logic(k, "ddg", point_arr, n, dim, 1);
    }
    else if (goal == 3){
        main_logic(k, "lnorm", point_arr, n, dim, 1);
    }
    else if (goal==4){
        main_logic(k, "jacobi", point_arr, n, dim, 1);
    }
    return NULL;
}


static PyObject* c_code_kmeans(PyObject *self, PyObject *args){
    PyObject *p_centroids, *p_observations, *item, *lst;
    int k, n, dim, i, j, k_item;
    double d_item;
    double **observations;
    double **final_centroids;
    int *ini_centroid_indices;
    struct Cluster** clusters;

    if (!PyArg_ParseTuple(args, "OOiii", &p_centroids, &p_observations, &k, &n, &dim)){
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
    assert(ini_centroid_indices!= NULL);

    for(i=0;i<k;i++){
        item = PyList_GetItem(p_centroids, i);
        k_item = PyLong_AsLong(item);
        ini_centroid_indices[i] = k_item;
    }

    clusters = (struct Cluster**) calloc(k, sizeof(struct Cluster*));
    assert(clusters != NULL);

    kmeans(observations, ini_centroid_indices, clusters, k, n, dim);

    free(observations);
    free(ini_centroid_indices);

    free(final_centroids);
    return NULL;
}

#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
        FUNC(METH_VARARGS, c_code, "first part of code"),
        FUNC(METH_VARARGS, c_code_kmeans, "kmeans part of code"),
        {NULL, NULL, 0, NULL}   /* sentinel */
};

 static struct PyModuleDef _moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        _methods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void){
    return PyModule_Create(&_moduledef);
}