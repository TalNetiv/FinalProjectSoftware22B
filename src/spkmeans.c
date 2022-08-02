#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int n, d; /*send to functions*/
void errorOccured();
static void emptyCentroids(int k, double** centroids);
static void emptyClusters(int k, double ***clusters);
static void addVectors(double *vectorA, double *vectorB, int d);
static void distribute(int n, int k, int d, double*** clusters, double **cPoints, double **centroids);
static void updateCentroids(int k, int d, double**centroids, double ***clusters);
static void insertToCluster(double ** cluster, double *vector, int d);
static int isVectorZero(double *vector, int d);
static int bestCluster(int point_index, int k, int d, double **cPoints, double **centroids);
static int convergentTrue(double ** curr, double ** prev, int d, int k);
double euclideanDistance(double* vectorA, double* vectorB, int d);
static void calculateCentroids(int cluster_index, int d, double ***clusters, double *centroid);
static void getPrevCentroids(int k, int d, double** centroids, double **getPrevCentroids);
static double*** createClusters(int k, int d, int n);
double ** matMultiply(double ** mat1, double ** mat2, int n);
double ** processPyObject(PyObject* data_points, int n, int d);
static PyObject *kMeans(PyObject *data_points, PyObject *initial_centroids, int n,int k,int d,double eps,int max_iter);
static PyObject *weightedAdjMat(PyObject *data_points, int n, int d);
static PyObject *diagDegMat(PyObject *data_points, int n, int d);
static PyObject *normalGraphLap(PyObject *data_points, int n, int d);
static PyObject *jacobian(PyObject *sym_max, int n);

void errorOccured() {
    printf("An Error Has Occured!");
}

static void emptyCentroids(int k, double **centroids) {
    int i;
    for (i = 0; i < k; i++) {
        free(centroids[i]);
    }
    free(centroids);
}

static void emptyClusters(int k, double ***clusters){
    int i;
    for (i = 0; i < k; i++) {
        free(clusters[i]);
    }
    free(clusters);
}

static void addVectors(double *vectorA, double *vectorB, int d){
    int i;
    for (i = 0; i < d; i++){
        vectorB[i] += vectorA[i];
    }
}

static void distribute(int n, int k, int d, double*** clusters, double **cPoints, double **centroids){
    int i;
    int index;
    for (i = 0; i < n; i++){
        index = bestCluster(i, k, d, cPoints, centroids);
        insertToCluster(clusters[index], cPoints[i], d);
    }
}

static void updateCentroids(int k, int d, double**centroids, double ***clusters){
    int i;
    for (i = 0; i < k; i++){
        calculateCentroids(i, d, clusters, centroids[i]);
    }
}

static void insertToCluster(double ** cluster, double *vector, int d){
    int i = 0;
    while (isVectorZero(cluster[i], d) == 0){
        i++;
    }
    cluster[i] = vector;
}

static int isVectorZero(double *vector, int d){
    int j;
    for (j = 0; j < d; j++){
        if (vector[j] != 0){
            return 0;
        }
    }
    return 1;
}

static int bestCluster(int point_index, int k, int d, double **cPoints, double **centroids){
    int m;
    int best_index = 0;
    double best_dist = pow(2,15);
    double curr_dist;
    for (m = 0; m < k; m++){
        curr_dist = euclideanDistance(cPoints[point_index], centroids[m], d);
        if (curr_dist < best_dist){
            best_dist = curr_dist;
            best_index = m;
        }
    }
    return best_index;
}

static int convergentTrue(double ** curr, double ** prev, int d, int k){
    int i;
    double eucdist;
    for (i=0; i<k; i++){
        eucdist = euclideanDistance(curr[i], prev[i], d); 
        if (eucdist > 0.001){
            return 0;
        }
    }
    return 1;
}

double euclideanDistance(double* vectorA, double* vectorB, int d){
    int i;
    double sum_squares = 0;
    for (i = 0 ; i < d ; i++){
        sum_squares += pow((vectorA[i]-vectorB[i]),2);
    }
    return sqrt(sum_squares);
}

static void calculateCentroids(int cluster_index, int d, double ***clusters, double *centroid){
    int i, p;
    int l = 0;
    int vector_index = 0;
    for (p = 0; p < d; p++){
        centroid[p] = 0;
    }
    while (isVectorZero(clusters[cluster_index][vector_index], d) != 1){
        addVectors(clusters[cluster_index][vector_index], centroid, d);
        l++;
        vector_index++;
    }
    for (i = 0; i < d; i++){
        if (l != 0){
            centroid[i] = centroid[i] / l;
    }   }
}

static void getPrevCentroids(int k, int d, double** centroids, double **prevCentroids){
    int i, j;
    for (i = 0; i < k; i++){
        for (j = 0; j < d; j++){  
            prevCentroids[i][j] = centroids[i][j];
        }
    }
}

static double*** createClusters(int k, int d, int n){
    int i, j;
    double*** new_clusters;
    new_clusters = (double***)calloc(k, sizeof(double **));
    if (new_clusters == NULL){
        errorOccured();
    }
    for (i = 0; i < k; i++){
        new_clusters[i] = (double **)calloc(n, sizeof(double *));
        if (new_clusters[i] == NULL){
            errorOccured();
        }
        for (j = 0; j < n; j++){
            new_clusters[i][j] = (double *)calloc(d, sizeof(double));
            if (new_clusters[i][j] == NULL) {
                errorOccured();
            }
        }
    }
    return new_clusters;
}

double ** matMultiply(double ** mat1, double ** mat2, int n) {
    int i, j, k;
    double sum = 0;
    double ** multResult = (double**)calloc(n, sizeof(double*));
    if (multResult == NULL) {
        errorOccured();
    }
    for (i = 0 ; i < n ; i ++) {
        multResult[i] = (double*)calloc(n, sizeof(double));
        if (multResult[i] == NULL) { errorOccured(); }
    }
    for (i = 0 ; i < n ; i ++) {
        for (j = 0; j < n; j ++){
            sum = 0;
        for (k = 0; k < n; k ++){
            sum = sum + mat1[i][k]*mat2[k][j];
        }
        multResult[i][j] = sum;
    }
}
    return multResult;
}

double ** processPyObject(PyObject *data_points, int n, int d) {
    int i, j;
    double **cPoints;
    cPoints = (double **)calloc(n, sizeof(double *)); /*process data points in seperate function*/
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

static PyObject* kMeans(PyObject *data_points, PyObject *initial_centroids, int n, int k, int d, double eps, int max_iter){
    int i, j, iter;
    double **centroids, **prevCentroids;
    double *** clusters;
    double ** cPoints = processPyObject(data_points, n, d);
    PyObject *res;
    centroids = (double **)calloc(k, sizeof(double *));
    if (centroids == NULL){
        errorOccured();
    }
    for (i = 0; i < k; i++){
        centroids[i] = (double *)calloc(d, sizeof(double));
        if (centroids[i] == NULL){
            errorOccured();
    }
        for (j = 0; j < d; j++){  
            centroids[i][j] = PyFloat_AsDouble(PyList_GetItem(initial_centroids,i*d +j));
        }
    }
    prevCentroids = (double **)calloc(k, sizeof(double *));
    if (prevCentroids == NULL){
        errorOccured();
    }
    for (i = 0; i < k; i++){
        prevCentroids[i] = (double *)calloc(d, sizeof(double));
        if(prevCentroids[i] == NULL){
            errorOccured();
        }
    }
    getPrevCentroids(k, d, centroids, prevCentroids);
    clusters = createClusters(k, d, n);
    distribute(n, k, d, clusters, cPoints, centroids);
    updateCentroids(k, d, centroids, clusters);
    iter = 0;
    while (convergentTrue(centroids, prevCentroids, d, k) == 0 && iter < max_iter){
        getPrevCentroids(k, d, centroids, prevCentroids);
        emptyClusters(k, clusters);
        clusters = createClusters(k, d, n);
        distribute(n, k, d, clusters, cPoints, centroids);
        updateCentroids(k, d, centroids, clusters);
        iter++;
    }
    res = PyList_New(0);
    PyObject *pylist;
    for (i = 0; i < k; ++i) {
        pylist = PyList_New(0);
        for (j = 0; j < d; ++j) {
            if (PyList_Append(pylist, PyFloat_FromDouble(centroids[i][j])) != 0)
            return NULL;
        }
        if (PyList_Append(res, pylist) != 0)
        return NULL;
    }

    emptyClusters(k, clusters);
    emptyCentroids(k, centroids);
    emptyCentroids(k, prevCentroids);
    emptyCentroids(n, cPoints);

    return res;
}
 
static PyObject* weightedAdjMat(PyObject *data_points, int n, int d){ /*wam */
    int i, j;
    double ** points = processPyObject(data_points, n, d);
    double** mat = (double**)calloc(n, sizeof(double*));
    if (mat == NULL) { errorOccured(); }
    for (i = 0 ; i < n ; i ++) {
        mat[i] = (double*)calloc(n, sizeof(double));
        if (mat[i] == NULL) { errorOccured(); }
    }
    for (i = 0 ; i < n ; i ++) {
        for (j = i ; j < (n) ; j++) {
            if (i == j) {
                mat[i][j] = 0;
            } else {
                double* v1 = points[i];
                double* v2 = points[j];
                double euc = euclideanDistance(v1, v2, d);
                double weight = exp((-euc)/2);
                mat[i][j] = weight;
                mat[j][i] = weight;
            }
        }
    }
    return mat;
}
 
static PyObject* diagDegMat(PyObject *data_points, int n, int d) { /* ddg */
    double sum;
    int i, j;
    double ** weights = weightedAdjMat(data_points, n, d);
    double** mat = (double**)calloc(n, sizeof(double*));
    if (mat == NULL) { errorOccured(); }
    for (i = 0 ; i < n ; i ++) {
        mat[i] = (double*)calloc(n, sizeof(double));
        if (mat[i] == NULL) { errorOccured(); }
    }
    for (i = 0 ; i < n ; i++) {
        sum = 0;
        for (j = 0; j < n ; j++) {
            sum += weights[i][j];
            }
        mat[i][i] = sum;
    }
    return mat;
}

static PyObject* normalGraphLap(PyObject *data_points, int n, int d) { /* lnorm */
    double ** multiplied;
    int i, j;
    double ** weights = weightedAdjMat(data_points, n, d);
    double ** diagmat = diagDegMat(data_points, n, d);
    double ** diagmatnew = (double**)calloc(n, sizeof(double*));
    if (diagmatnew == NULL) { errorOccured(); }
    for (i = 0 ; i < n ; i ++) {
        diagmatnew[i] = (double*)calloc(n, sizeof(double));
        if (diagmatnew[i] == NULL) { errorOccured(); }
    }
    for (i = 0 ; i < n ; i ++) {
        diagmatnew[i][i] = (1/sqrt(diagmat[i][i]));
    }
    double ** mat = (double**)calloc(n, sizeof(double*));
    if (mat == NULL) { errorOccured(); }
    for (i = 0 ; i < n ; i ++) {
        mat[i] = (double*)calloc(n, sizeof(double));
        if (mat[i] == NULL) { errorOccured(); }
    }
    multiplied = matMultiply(matMultiply(diagmatnew, weights, n), diagmatnew, n);
    for (i = 0 ; i < n ; i++) {
        for (j = 0 ; j < n ; j++) {
            if (i == j) {
                mat[i][j] = 1 - multiplied[i][j];
            } else {
            mat[i][j] = -multiplied[i][j];
            }
        }
    }
    return mat;
}

static PyObject *jacobian(PyObject *sym_max, int n) {

}

