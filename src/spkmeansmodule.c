#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.c"

/* args are the arguments passed from the python program*/
static PyObject* spk(PyObject *self, PyObject *args){
    PyObject *data_points, *initial_centroids; /* both in the format of long doubles list */
    int n, k, d, max_iter;
    double eps;

    if (!PyArg_ParseTuple(args, "OOiiiid", &data_points, &initial_centroids, &n, &k, &d, &max_iter, &eps)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    return kMeans(data_points, initial_centroids, n, k, d, eps, max_iter);
}

static PyObject* wam(PyObject *self, PyObject *args){
    PyObject *data_points; /* both in the format of long doubles list */
    int n, d;

    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &d)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    return weightedAdjMat(data_points, n, d);
}

static PyObject* ddg(PyObject *self, PyObject *args){
    PyObject *data_points; /* both in the format of long doubles list */
    int n, d;

    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &d)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    return diagDegMat(data_points, n, d);
}

static PyObject* lnorm(PyObject *self, PyObject *args){
    PyObject *data_points; /* both in the format of long doubles list */
    int n, d;

    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &d)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    return normalGraphLap(data_points, n, d);
}

static PyObject* jacobi(PyObject *self, PyObject *args){
    PyObject *sym_mat; /* both in the format of long doubles list */
    int n;

    if (!PyArg_ParseTuple(args, "Oi", &sym_mat, &n)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    return jacobian(sym_mat, n);
}

static PyMethodDef capiMethods[] = {
    {"spk", /* name to be used */
        (PyCFunction) spk, /* the C function that implements and returns statis PyObject* */
        METH_VARARGS,
        PyDoc_STR("The function receives list of initial centroids and data points, and returns final clustered centroids")
    },
    {"wam", /* name to be used */
    (PyCFunction) wam, /* the C function that implements and returns statis PyObject* */
    METH_VARARGS,
    PyDoc_STR("Calculate and output the Weighted Adjacency Matrix")
    },
    {"ddg", /* name to be used */
    (PyCFunction) ddg, /* the C function that implements and returns statis PyObject* */
    METH_VARARGS,
    PyDoc_STR("Calculate and output the Diagonal Degree Matrix")
    },
    {"lnorm", /* name to be used */
    (PyCFunction) lnorm, /* the C function that implements and returns statis PyObject* */
    METH_VARARGS,
    PyDoc_STR("Calculate and output the Normalized Graph Laplacian")
    },
    {"jacobi", /* name to be used */
    (PyCFunction) jacobi, /* the C function that implements and returns statis PyObject* */
    METH_VARARGS,
    PyDoc_STR("Calculate and output the eigenvalues and eigenvectors")
    },
    {NULL, NULL, 0, NULL}
};
/* indicates the module using the above definitions */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmodule", /* name of module */
    NULL, /* empty module documentation */
    -1,
    capiMethods /* array from above */
};

PyMODINIT_FUNC
PyInit_spkmodule(void) { /* creates the module, this should be the only non-static item in this file */
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
