#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

double ** PyToC(PyObject *data_points, int n, int d) {
    int i, j;
    double **cPoints;
    cPoints = (double **)calloc(n, sizeof(double *));
    if (cPoints == NULL){
        errorOccured();
    }
    for (i = 0; i < n; i++){
        cPoints[i] = (double *)calloc(d, sizeof(double));
        if (cPoints[i] == NULL){
            errorOccured();
    }
    }
    for (i = 0; i < n; i++){
        for (j = 0; j < d; j++){
            cPoints[i][j]=PyFloat_AsDouble(PyList_GetItem(data_points, i*d + j));
        }
    } 
    return cPoints;
}

static PyObject* CToPy(double** mat, int n, int d) {
    int i, j;
    PyObject * res;
    res = PyList_New(0);
    PyObject *pylist;
    for (i = 0; i < n; i++) {
        pylist = PyList_New(0);
        for (j = 0; j < d; j++) {
            if (PyList_Append(pylist, PyFloat_FromDouble(mat[i][j])) != 0)
            return NULL;
        }
        if (PyList_Append(res, pylist) != 0)
        return NULL;
    }
    return res;
}

/* args are the arguments passed from the python program*/
static PyObject* spk(PyObject *self, PyObject *args){
    PyObject *data_points, *initial_centroids; /* both in the format of long doubles list */
    int n, k, d;
    double** points, ** centroids;

    if (!PyArg_ParseTuple(args, "OOiii", &data_points, &initial_centroids, &n, &k, &d)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    points = PyToC(data_points, n, d);
    centroids = PyToC(initial_centroids, k, d);
    return CToPy(kMeans(points, centroids, n, k, d), k, d);
}

static PyObject* wam(PyObject *self, PyObject *args){
    PyObject *data_points; /* both in the format of long doubles list */
    int n, d;
    double** points;

    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &d)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }

    
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    points = PyToC(data_points, n, d);
    return CToPy(weightedAdjMat(points, n, d), n, n);
}

static PyObject* ddg(PyObject *self, PyObject *args){
    PyObject *data_points; /* both in the format of long doubles list */
    int n, d;
    double** points;

    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &d)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    points = PyToC(data_points, n, d);
    return CToPy(diagDegMat(points, n, d), n, n);
}

static PyObject* lnorm(PyObject *self, PyObject *args){
    PyObject *data_points; /* both in the format of long doubles list */
    int n, d;
    double** points;

    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &d)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    points = PyToC(data_points, n, d);
    return CToPy(normalGraphLap(points, n, d), n, n);
}

static PyObject* jacobi(PyObject *self, PyObject *args){
    PyObject *sym_mat; /* both in the format of long doubles list */
    int n;
    double** entries;

    if (!PyArg_ParseTuple(args, "Oi", &sym_mat, &n)){
        return NULL; /*  NULL implies an error occured because it's not allowed for PyObject* to be NULL*/
    }
    /* build the answer ("d" = convert a C double to a python floating point number) back into a python object */
    entries = PyToC(sym_mat, n, n);
    return CToPy(jacobian(entries, n), n+1, n);
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
