#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int n, d; /*send to functions*/
static void errorOccured();
static double ** initializeMat(int n, int d);
static double ** readFromFile(char* fileName);
static void printMat(double** mat, int n);
static void emptyCentroids(int k, double** centroids);
static void emptyClusters(int k, double ***clusters);
static void addVectors(double *vectorA, double *vectorB, int d);
static void distribute(int n, int k, int d, double*** clusters, double **cPoints, double **centroids);
static void updateCentroids(int k, int d, double**centroids, double ***clusters);
static void insertToCluster(double ** cluster, double *vector, int d);
static int isVectorZero(double *vector, int d);
static int bestCluster(int point_index, int k, int d, double **cPoints, double **centroids);
static int convergentTrue(double ** curr, double ** prev, int d, int k);
static double euclideanDistance(double* vectorA, double* vectorB, int d);
static void calculateCentroids(int cluster_index, int d, double ***clusters, double *centroid);
static void getPrevCentroids(int k, int d, double** centroids, double **getPrevCentroids);
static double*** createClusters(int k, int d, int n);
static double** matMultiply(double ** mat1, double ** mat2, int n);
static int* findGreatestValue(double** mat, int n);
static int isMatDiag(double** mat, int n);
static double off(double** mat, int n);
static double* obtainVariables(double** mat, int n, int row, int col);
static double** createRotMat(int n, double c, double s, int row, int col);
static double** kMeans(double** points, double** initial_centroids, int n, int k, int d);
static  double** weightedAdjMat(double** data_points, int n, int d);
static double** diagDegMat(double** data_points, int n, int d);
static double** normalGraphLap(double** data_points, int n, int d);
/* static double** jacobian(double** sym_max, int n); */

static void errorOccured() {
    printf("An Error Has Occured!");
    exit(1);
}

static double ** initializeMat(int n, int d){
    int i;
    double ** mat = (double**)calloc(n, sizeof(double*));
    if (mat == NULL) { errorOccured(); }
    for (i = 0 ; i < n ; i ++) {
        mat[i] = (double*)calloc(d, sizeof(double));
        if (mat[i] == NULL) { errorOccured(); }
    }
    return mat;
}

static double ** readFromFile(char* fileName){
    double ** points;
    int i, j;
    double cord;
    char psik, c;
    FILE *ifp;
    ifp = fopen(fileName, "r");
    if (ifp == NULL) {
        errorOccured();
    }
    printf("passed first error barrier");
    c = fgetc(ifp);
    while (c !='\n'){
        if (c == ','){
            d++;
        }
        c = fgetc(ifp);
    }
    d++;
    while (c != EOF){
        if (c == '\n'){
        n++;
        }
        c = fgetc(ifp);
    }
    rewind(ifp);
    printf("go initialize points in read from file");
    points = initializeMat(n, d);
    for (i = 0; i < n; i++){
        for (j = 0; j < d; j++){
            fscanf(ifp, "%lf%c",&cord, &psik);
            points[i][j] = cord;
        }
    }
    fclose(ifp);
    return points;
}

static void printMat(double** mat, int n){
    int rows, columns;
 	for(rows = 0; rows < n; rows++)
  	{
        printf("row num: %d  " , rows); /* delete on submission !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  		for(columns = 0; columns < n; columns++)
  		{
  			if (columns == n - 1){
                printf("%.4f \n" , mat[rows][columns]);
            } else{
                printf("%.4f,", mat[rows][columns]);
            }
		}
  	}  	
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
    double ** multResult = initializeMat(n, n);
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

static int* findGreatestValue(double** mat, int n) {
    int i, j, row, col;
    double maxval = 0;
    for (i = 1; i < n; i++) {
        for (j = i; j < n; j++) {
            if (abs(mat[i][j]) > maxval) {
                maxval = abs(mat[i][j]);
                row = i;
                col = j;
            }
        }
    }
    int* entry = {row, col};
    return entry;
}

static int isMatDiag(double** mat, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j && mat[i][j] != 0) {
                return 0;
            }
        }
    }
    return 1;
}

static double off(double** mat, int n) {
    int i, j;
    double sum;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                sum += pow(mat[i][j], 2);
            }
        }
    }
    return sum;
}

static double* obtainVariables(double** mat, int n, int row, int col) {
    double theta, c, t;
    double* variables = (double*)calloc(3, sizeof(double));
    int sign;
    if (variables == NULL) {errorOccured();}
    theta = (mat[row][row]-mat[col][col])/(2*mat[row][col]);
    variables[0] = theta;
    if (sign < 0) {sign = -1;}
    else {sign = 1;}
    t = sign/(abs(theta)+sqrt(pow(theta, 2)+1));
    variables[1] = t;
    c = 1/(sqrt(pow(t, 2)+1));
    variables[2] = c;
    return variables;
}

static double** createRotMat(int n, double c, double s, int row, int col) {
    int i, j;
    double** P = initializeMat(n, n);
    for (i = 0; i < n; i++) {
        P[i][i] = 1;
    }
    P[row][row] = c;
    P[col][col] = c;
    P[row][col] = s;
    P[col][row] = -s;
    return P;
}

static double** kMeans(double** points, double** initial_centroids, int n, int k, int d){
    int iter;
    double **centroids, **prevCentroids;
    double *** clusters;
    centroids = initializeMat(k, d);
    centroids = initial_centroids;
    prevCentroids = initializeMat(k, d);
    getPrevCentroids(k, d, centroids, prevCentroids);
    clusters = createClusters(k, d, n);
    distribute(n, k, d, clusters, points, centroids);
    updateCentroids(k, d, centroids, clusters);
    iter = 0;
    while (convergentTrue(centroids, prevCentroids, d, k) == 0 && iter < 300){
        getPrevCentroids(k, d, centroids, prevCentroids);
        emptyClusters(k, clusters);
        clusters = createClusters(k, d, n);
        distribute(n, k, d, clusters, points, centroids);
        updateCentroids(k, d, centroids, clusters);
        iter++;
    }
    emptyClusters(k, clusters);
    emptyCentroids(k, centroids);
    emptyCentroids(k, prevCentroids);
    emptyCentroids(n, points);
    return centroids;
}

 
static double** weightedAdjMat(double** points, int n, int d){ /*wam */
    int i, j;
    double** mat = initializeMat(n, n);
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
 
static double** diagDegMat(double** points, int n, int d) { /* ddg */
    double sum;
    int i, j;
    double ** weights = weightedAdjMat(points, n, d);
    double** mat = initializeMat(n, n);
    for (i = 0 ; i < n ; i++) {
        sum = 0;
        for (j = 0; j < n ; j++) {
            sum += weights[i][j];
            }
        mat[i][i] = sum;
    }
    return mat;
}

static double** normalGraphLap(double** points, int n, int d) { /* lnorm */
    double ** multiplied;
    int i, j;
    double ** weights = weightedAdjMat(points, n, d);
    double ** diagmat = diagDegMat(points, n, d);
    double ** diagmatnew = initializeMat(n, n);
    double ** mat = initializeMat(n, n);
    for (i = 0 ; i < n ; i ++) {
        diagmatnew[i][i] = (1/sqrt(diagmat[i][i]));
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

static double** jacobian(double** A, int n) {
    int i, j, iter, piv_row, piv_col;
    double phi, c, s;
    double* variables;
    double** Atag, **pivot, **P, **Ptranspose;
    double eps = 1.0*pow(10, -5);
    Atag = A;
    do { 
        iter++;
        A = Atag;
        pivot = findGreatestValue(A,n);
        piv_row = pivot[0];
        piv_col = pivot[1];
        variables = obtainVariables(A, n, piv_row, piv_col);
        phi = 1/(2*tan(1/variables[0]));
        c = variables[2];
        s = c*variables[1];
        P = createRotMat(n, c, s, piv_row, piv_col);
        Ptranspose = initializeMat(n, n);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                Ptranspose[j][i] = P[i][j];
            }
        }
        Atag = matMultiply(matMultiply(Ptranspose, A, n), P, n);
        if (isMatDiag(Atag, n) == 1) {
            break;
        }
    } while ((off(A, n) - off(Atag, n) > eps) && iter < 100);

}

int main(int argc, char *argv[]){

    double **points, **mat;
    char *goal = argv[1];
    char *fileName = argv[2];

    if (argc != 3){ 
        printf("Invalid Input!Main");
        exit(1);
    }
    if ( (strcmp(goal, "wam") != 0) && (strcmp(goal, "ddg") != 0) && (strcmp(goal, "lnorm") != 0) && (strcmp(goal, "jacobi") != 0) ){  /*check goal validity*/ 
        printf("Invalid Input!Main");
        exit(1);
    }

    printf("go read from file");
    points = readFromFile(fileName);
    printf("done reading from file\n");

 /*   if (strcmp(goal, "jacobi") == 0){  ###add free mat func
        res = initializeMat(n+1, n);
        res = jacobi(points);
        print2Darray(res, n, n);
        free2D(points);
        free2D(res);
        return 0;
    } */
    mat = initializeMat(n, n);
    printf("ken");
    if (strcmp(goal, "wam") == 0){
        mat = weightedAdjMat(points, n, d);
    }
    if (strcmp(goal, "ddg") == 0){
        mat = diagDegMat(points, n, d);
    }
    if (strcmp(goal, "lnorm") == 0){
        mat = normalGraphLap(points, n, d);
    }
    if (strcmp(goal, "charta") == 0){
        mat = kMeans(mat, mat, 1,1,1);
    }
    /* add free to matrices here */
    printMat(mat, n);
    return 0; 
}