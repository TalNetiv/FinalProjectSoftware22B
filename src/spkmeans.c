#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int n, d;
static void errorOccured();
static void invalidInput();
static double **initializeMat(int n, int d);
static double **readFromFile(char *fileName);
static void printMat(double **mat, int n, int d);
static void emptyMatrix(int k, double **matrix);
static void emptyClusters(int k, double ***clusters);
static void addVectors(double *vectorA, double *vectorB, int d);
static void distribute(int n, int k, int d, double ***clusters, double **cPoints, double **centroids);
static void updateCentroids(int k, int d, double **centroids, double ***clusters);
static void insertToCluster(double **cluster, double *vector, int d);
static int isVectorAllocated(double *vector, int d);
static int bestCluster(int point_index, int k, int d, double **cPoints, double **centroids);
static int convergentTrue(double **curr, double **prev, int d, int k);
static double euclideanDistance(double *vectorA, double *vectorB, int d);
static void calculateCentroids(int cluster_index, int d, double ***clusters, double *centroid);
static void getPrevCentroids(int k, int d, double **centroids, double **getPrevCentroids);
static double ***createClusters(int k, int d, int n);
static double **matMultiply(double **mat1, double **mat2, int n);
static int *findGreatestValue(double **mat, int n);
static int isMatDiag(double **mat, int n);
static double off(double **mat, int n);
static double *obtainVariables(double **mat, int row, int col);
static double **createRotMat(int n, double c, double s, int row, int col);
static double **kMeans(double **points, double **initial_centroids, int n, int k, int d);
static double **weightedAdjMat(double **points, int n, int d);
static double **diagDegMat(double **points, int n, int d);
static double **normalGraphLap(double **points, int n, int d);
static double **jacobian(double **A, int n);

static void errorOccured()
{
    printf("An Error Has Occured");
    exit(1);
}

static void invalidInput()
{
    printf("Invalid Input!");
    exit(1);
}

/* A function that allocates space for an nxd matrix */
static double **initializeMat(int n, int d)
{
    int i;
    double **mat = (double **)calloc(n, sizeof(double *));
    if (mat == NULL)
    {
        errorOccured();
    }
    for (i = 0; i < n; i++)
    {
        mat[i] = (double *)calloc(d, sizeof(double));
        if (mat[i] == NULL)
        {
            errorOccured();
        }
    }
    return mat;
}

static double **readFromFile(char *fileName)
{
    double **points;
    int i, j;
    double cord;
    char psik, c;
    FILE *ifp;
    ifp = fopen(fileName, "r");
    if (ifp == NULL)
    {
        errorOccured();
    }
    c = fgetc(ifp);
    while (c != '\n')
    {
        if (c == ',')
        {
            d++; /* d is the amount of entries each data point has */
        }
        c = fgetc(ifp);
    }
    d++;
    while (c != EOF)
    {
        if (c == '\n')
        {
            n++; /* n is the amount of data points */
        }
        c = fgetc(ifp);
    }
    rewind(ifp);
    points = initializeMat(n, d);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            fscanf(ifp, "%lf%c", &cord, &psik);
            points[i][j] = cord;
        }
    }
    fclose(ifp);
    return points;
}

static void printMat(double **mat, int n, int d)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            if (j == d - 1)
            {
                printf("%.4f\n", mat[i][j]);
            }
            else
            {
                printf("%.4f,", mat[i][j]);
            }
        }
    }
}

/* free allocated space */
static void emptyMatrix(int k, double **matrix)
{
    int i;
    for (i = 0; i < k; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

/* empty all clusters in the kmeans algorithm before redistributing into them */
static void emptyClusters(int k, double ***clusters)
{
    int i;
    for (i = 0; i < k; i++)
    {
        free(clusters[i]);
    }
    free(clusters);
}

static void addVectors(double *vectorA, double *vectorB, int d)
{
    int i;
    for (i = 0; i < d; i++)
    {
        vectorB[i] += vectorA[i];
    }
}

/* for every point, find the cluster represented by the closest centroid and insert the point in there */
static void distribute(int n, int k, int d, double ***clusters, double **cPoints, double **centroids)
{
    int i;
    int index;
    for (i = 0; i < n; i++)
    {
        index = bestCluster(i, k, d, cPoints, centroids);
        insertToCluster(clusters[index], cPoints[i], d);
    }
}

/* after distributing the points, update each cluster's centroid */
static void updateCentroids(int k, int d, double **centroids, double ***clusters)
{
    int i;
    for (i = 0; i < k; i++)
    {
        calculateCentroids(i, d, clusters, centroids[i]);
    }
}

static void insertToCluster(double **cluster, double *vector, int d)
{
    int i = 0;
    while (isVectorAllocated(cluster[i], d) == 0)
    {
        i++;
    }
    cluster[i] = vector;
}

/* in every iteration we don't know how many points would be allocated to each cluster, thus we initialize all clusters to
 the size of nXd, and this function checks if the current vector is an actual point or merely a space allocated */
static int isVectorAllocated(double *vector, int d)
{
    int j;
    for (j = 0; j < d; j++)
    {
        if (vector[j] != -210496)
        {
            return 0;
        }
    }
    return 1;
}

/* function to be used when distributing points to clusters */
static int bestCluster(int point_index, int k, int d, double **cPoints, double **centroids)
{
    int m;
    int best_index = 0;
    double best_dist = pow(2, 20);
    double curr_dist;
    for (m = 0; m < k; m++)
    {
        curr_dist = euclideanDistance(cPoints[point_index], centroids[m], d);
        if (curr_dist < best_dist)
        {
            best_dist = curr_dist;
            best_index = m;
        }
    }
    return best_index;
}

/* this function checks if the kmeans convergence condition has been reached */
static int convergentTrue(double **curr, double **prev, int d, int k)
{
    int i;
    double eucdist;
    for (i = 0; i < k; i++)
    {
        eucdist = euclideanDistance(curr[i], prev[i], d);
        if (eucdist > 0)
        {
            return 0;
        }
    }
    return 1;
}

double euclideanDistance(double *vectorA, double *vectorB, int d)
{
    int i;
    double sum_squares = 0;
    for (i = 0; i < d; i++)
    {
        sum_squares += pow((vectorA[i] - vectorB[i]), 2);
    }
    return sqrt(sum_squares);
}

/* this function computes the new centroid for each cluster after distributing points to the clusters */
static void calculateCentroids(int cluster_index, int d, double ***clusters, double *centroid)
{
    int i, p;
    int l = 0;
    int vector_index = 0;
    for (p = 0; p < d; p++)
    {
        centroid[p] = 0;
    }
    while (isVectorAllocated(clusters[cluster_index][vector_index], d) != 1)
    {
        addVectors(clusters[cluster_index][vector_index], centroid, d);
        l++;
        vector_index++;
    }
    for (i = 0; i < d; i++)
    {
        if (l != 0)
        {
            centroid[i] = centroid[i] / l;
        }
    }
}

static void getPrevCentroids(int k, int d, double **centroids, double **prevCentroids)
{
    int i, j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            prevCentroids[i][j] = centroids[i][j];
        }
    }
}

/* initialize k clusters with junk values */
static double ***createClusters(int k, int d, int n)
{
    int i, j, l;
    double ***new_clusters;
    new_clusters = (double ***)calloc(k, sizeof(double **));
    if (new_clusters == NULL)
    {
        errorOccured();
    }
    for (i = 0; i < k; i++)
    {
        new_clusters[i] = (double **)calloc(n, sizeof(double *));
        if (new_clusters[i] == NULL)
        {
            errorOccured();
        }
        for (j = 0; j < n; j++)
        {
            new_clusters[i][j] = (double *)calloc(d, sizeof(double));
            if (new_clusters[i][j] == NULL)
            {
                errorOccured();
            }
            for (l = 0; l < d; l++)
            {
                new_clusters[i][j][l] = -210496;
            }
        }
    }
    return new_clusters;
}

/* multiplies nXn matrices */
double **matMultiply(double **mat1, double **mat2, int n)
{
    int i, j, k;
    double sum = 0;
    double **multResult = initializeMat(n, n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            sum = 0;
            for (k = 0; k < n; k++)
            {
                sum = sum + mat1[i][k] * mat2[k][j];
            }
            multResult[i][j] = sum;
        }
    }
    return multResult;
}

/* This function returns the row and the column of the greastest entry (in absolute value)
in the off-diagonal part of a matrix */
static int *findGreatestValue(double **mat, int n)
{
    int i, j;
    double maxval = 0;
    int *entry = (int *)calloc(2, sizeof(int));
    if (entry == NULL)
    {
        errorOccured();
    }
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            if (fabs(mat[i][j]) > maxval)
            {
                maxval = fabs(mat[i][j]);
                entry[0] = i;
                entry[1] = j;
            }
        }
    }
    if (maxval == 0)
    {
        entry[0] = 0;
        entry[1] = 1;
    }
    return entry;
}

static int isMatDiag(double **mat, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i != j && mat[i][j] != 0)
            {
                return 0;
            }
        }
    }
    return 1;
}

/* used in jacobi algorithm convergence condition check */
static double off(double **mat, int n)
{
    int i, j;
    double sum = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                sum += pow(mat[i][j], 2);
            }
        }
    }
    return sum;
}

/* Using the definitions of theta, t and c in the Obtain section
to find them using the greatest entry of the matrix */
static double *obtainVariables(double **mat, int row, int col)
{
    double theta, c, t;
    int sign;
    double *variables = (double *)calloc(3, sizeof(double));
    if (variables == NULL)
    {
        errorOccured();
    }
    if (mat[row][col] == 0)
    {
        variables[0] = strtod("Inf", NULL);
        variables[1] = 0;
        variables[2] = 1;
        return variables;
    }
    theta = (mat[col][col] - mat[row][row]) / (2 * mat[row][col]);
    variables[0] = theta;
    if (theta < 0)
    {
        sign = -1;
    }
    else
    {
        sign = 1;
    }
    t = sign / (fabs(theta) + pow((pow(theta, 2) + 1), 0.5));
    variables[1] = t;
    c = pow(pow(t, 2) + 1, -0.5);
    variables[2] = c;
    return variables;
}

/* creates the P rotation matrix for jacobi algorithm */
static double **createRotMat(int n, double c, double s, int row, int col)
{
    int i;
    double **P = initializeMat(n, n);
    for (i = 0; i < n; i++)
    {
        P[i][i] = 1;
    }
    P[row][row] = c;
    P[col][col] = c;
    P[row][col] = s;
    P[col][row] = s == 0 ? 0 : -s;
    return P;
}

static double **kMeans(double **points, double **initial_centroids, int n, int k, int d)
{
    int iter, i, j;
    double **centroids, **prevCentroids;
    double ***clusters;
    centroids = initializeMat(k, d);
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            centroids[i][j] = initial_centroids[i][j];
        }
    }
    prevCentroids = initializeMat(k, d);
    getPrevCentroids(k, d, centroids, prevCentroids);
    clusters = createClusters(k, d, n);
    distribute(n, k, d, clusters, points, centroids);
    updateCentroids(k, d, centroids, clusters);
    iter = 0;
    while (convergentTrue(centroids, prevCentroids, d, k) == 0 && iter < 300)
    {
        getPrevCentroids(k, d, centroids, prevCentroids);
        emptyClusters(k, clusters);
        clusters = createClusters(k, d, n);
        distribute(n, k, d, clusters, points, centroids);
        updateCentroids(k, d, centroids, clusters);
        iter++;
    }
    emptyClusters(k, clusters);
    emptyMatrix(k, prevCentroids);
    emptyMatrix(n, points);
    return centroids;
}

static double **weightedAdjMat(double **points, int n, int d)
{ /* wam */
    int i, j;
    double **mat = initializeMat(n, n);
    for (i = 0; i < n; i++)
    {
        for (j = i; j < (n); j++)
        {
            if (i == j)
            {
                mat[i][j] = 0;
            }
            else
            {
                double *v1 = points[i];
                double *v2 = points[j];
                double euc = euclideanDistance(v1, v2, d);
                double weight = exp((-euc) / 2);
                mat[i][j] = weight;
                mat[j][i] = weight;
            }
        }
    }
    return mat;
}

static double **diagDegMat(double **points, int n, int d)
{ /* ddg */
    double sum;
    int i, j;
    double **weights = weightedAdjMat(points, n, d);
    double **mat = initializeMat(n, n);
    for (i = 0; i < n; i++)
    {
        sum = 0;
        for (j = 0; j < n; j++)
        {
            sum += weights[i][j];
        }
        mat[i][i] = sum;
    }
    emptyMatrix(n, weights);
    return mat;
}

static double **normalGraphLap(double **points, int n, int d)
{ /* lnorm */
    double **multiplied;
    int i, j;
    double **weights = weightedAdjMat(points, n, d);
    double **diagmat = diagDegMat(points, n, d);
    double **diagmatnew = initializeMat(n, n);
    double **mat = initializeMat(n, n);
    for (i = 0; i < n; i++)
    {
        diagmatnew[i][i] = (1 / sqrt(diagmat[i][i]));
    }
    multiplied = matMultiply(matMultiply(diagmatnew, weights, n), diagmatnew, n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                mat[i][j] = 1 - multiplied[i][j];
            }
            else
            {
                mat[i][j] = -multiplied[i][j];
            }
        }
    }
    emptyMatrix(n, weights);
    emptyMatrix(n, diagmat);
    emptyMatrix(n, diagmatnew);
    return mat;
}

static double **jacobian(double **A, int n)
{ /* TODO: add special treatment for singleton matrix */
    int i, j, piv_row, piv_col;
    double c, s;
    double *variables;
    double **Atag, **P, **Ptranspose, **V, **mat;
    int *pivot;
    int iter = 0;
    Atag = A;
    while (iter <= 100)
    {
        iter++;
        A = Atag;
        pivot = findGreatestValue(A, n);
        piv_row = pivot[0];
        piv_col = pivot[1];
        variables = obtainVariables(A, piv_row, piv_col);
        c = variables[2];
        s = c * variables[1];
        P = createRotMat(n, c, s, piv_row, piv_col);
        if (iter == 1)
        {
            V = P;
        }
        else
        {
            V = matMultiply(V, P, n);
        }
        Ptranspose = initializeMat(n, n);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Ptranspose[i][j] = P[j][i];
            }
        }
        Atag = matMultiply(matMultiply(Ptranspose, A, n), P, n);
        if (isMatDiag(Atag, n) || (off(A, n) - off(Atag, n)) < 0.00001)
        {
            break;
        }
    }
    mat = initializeMat(n + 1, n);
    for (j = 0; j < n; j++)
    {
        mat[0][j] = Atag[j][j];
    }
    for (i = 1; i < n + 1; i++)
    {
        for (j = 0; j < n; j++)
        {
            mat[i][j] = V[i - 1][j];
        }
    }
    emptyMatrix(n, P);
    emptyMatrix(n, Ptranspose);
    emptyMatrix(n, A);
    emptyMatrix(n, Atag);
    return mat;
}

int main(int argc, char *argv[])
{
    double **points, **mat;
    char *goal, *fileName;
    if (argc != 3)
    {
        invalidInput();
    }
    goal = argv[1];
    fileName = argv[2];
    if ((strcmp(goal, "wam") != 0) && (strcmp(goal, "ddg") != 0) && (strcmp(goal, "lnorm") != 0) && (strcmp(goal, "jacobi") != 0))
    {
        invalidInput();
    }
    points = readFromFile(fileName);
    if (strcmp(goal, "jacobi") == 0)
    {
        mat = initializeMat(n + 1, n);
        mat = jacobian(points, n);
        printMat(mat, n + 1, n);
        emptyMatrix(n + 1, mat);
        return 0;
    }
    mat = initializeMat(n, n);
    if (strcmp(goal, "wam") == 0)
    {
        mat = weightedAdjMat(points, n, d);
    }
    if (strcmp(goal, "ddg") == 0)
    {
        mat = diagDegMat(points, n, d);
    }
    if (strcmp(goal, "lnorm") == 0)
    {
        mat = normalGraphLap(points, n, d);
    }
    /* this is to avoid compilation error since kMeans() is only called from module, and then not being used locally */
    if (strcmp(goal, "kmeansalgorithm") == 0)
    {
        mat = kMeans(points, points, 1, 1, 1);
    }
    printMat(mat, n, n);
    emptyMatrix(n, mat);
    return 0;
}
