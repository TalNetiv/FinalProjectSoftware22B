#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int dpamount, dplen, i, j, k;
double point;
double** dps;
char c, psik;
void InvalidInput();
void ErrorOccured();
void DataPointsData(char* name);
double euclideanDistance(double* vectorA, double* vectorB);
double** WeightedAdjMat();
double** DiagDegMat(double** weights);
double ** MatMultiply(double ** mat1, double ** mat2);
double** NormalGraphLap();

void ErrorOccured() {
    printf("An Error Has Occured!");
}

double euclideanDistance(double* vectorA, double* vectorB){
    int i;
    double sum_squares = 0;
    for (i = 0 ; i < dplen ; i++){
        sum_squares += pow((vectorA[i]-vectorB[i]),2);
    }
    return sqrt(sum_squares);
}

double** WeightedAdjMat() { /* wam */
    double** mat = (double**)calloc(dpamount, sizeof(double*));
    if (mat == NULL) { ErrorOccured(); }
    for (i = 0 ; i < dpamount ; i ++) {
        mat[i] = (double*)calloc(dpamount, sizeof(double));
        if (mat[i] == NULL) { ErrorOccured(); }
    }
    for (i = 0 ; i < dpamount ; i ++) {
        for (j = i ; j < (dpamount) ; j++) {
            if (i == j) {
                mat[i][j] = 0;
            } else {
                double* v1 = dps[i];
                double* v2 = dps[j];
                double euc = euclideanDistance(v1, v2);
                double weight = exp((-euc)/2);
                mat[i][j] = weight;
                mat[j][i] = weight;
            }
        }
    }
    return mat;
}

double** DiagDegMat(double** weights) { /* ddg */
    double sum;
    double** mat = (double**)calloc(dpamount, sizeof(double*));
    if (mat == NULL) { ErrorOccured(); }
    for (i = 0 ; i < dpamount ; i ++) {
        mat[i] = (double*)calloc(dpamount, sizeof(double));
        if (mat[i] == NULL) { ErrorOccured(); }
    }
    for (i = 0 ; i < dpamount ; i++) {
        sum = 0;
        for (k = 0; k < dpamount ; k++) {
            sum += weights[i][k];
            }
        mat[i][i] = sum;
    }
    return mat;
}

double ** MatMultiply(double ** mat1, double ** mat2) {

}

double ** NormalGraphLap(double ** weights, double ** diagmat) {
    double ** multiplied;
    double ** diagmatnew = (double**)calloc(dpamount, sizeof(double*));
    if (diagmatnew == NULL) { ErrorOccured(); }
    for (i = 0 ; i < dpamount ; i ++) {
        diagmatnew[i] = (double*)calloc(dpamount, sizeof(double));
        if (diagmatnew[i] == NULL) { ErrorOccured(); }
    }
    for (i = 0 ; i < dpamount ; i ++) {
        diagmatnew[i][i] = (1/sqrt(diagmat[i][i]));
    }
    double ** mat = (double**)calloc(dpamount, sizeof(double*));
    if (mat == NULL) { ErrorOccured(); }
    for (i = 0 ; i < dpamount ; i ++) {
        mat[i] = (double*)calloc(dpamount, sizeof(double));
        if (mat[i] == NULL) { ErrorOccured(); }
    }
    multiplied = MatMultiply(MatMultiply(diagmatnew, weights), diagmatnew);
    for (i = 0 ; i < dpamount ; i++) {
        for (j = 0 ; j < dpamount ; j++) {
            if (i==j) {
                mat[i][j] = 1 - multiplied[i][j];
            } else {
            mat[i][j] = -multiplied[i][j];
            }
        }
    }
    return mat;
}