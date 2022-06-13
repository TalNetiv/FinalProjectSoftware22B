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
double** NormalGraphLap();

void InvalidInput() {
    printf("Invalid Input!");
}

void ErrorOccured() {
    printf("An Error Has Occured!");
}

void DataPointsData(char* name) {
    FILE * input = fopen(name, "r");
    if (input == NULL) {
        ErrorOccured();
    }
    dplen = 0, dpamount = 0;
    c = fgetc(input);
    while (c !='\n'){
        if (c == ','){
            dplen++;
        }
        c = fgetc(input);
    }
    dplen++;
    while (c != EOF){
        if (c == '\n'){
            dpamount++;
        }
        c = fgetc(input);
    }
    rewind(input);
    dps = (double**)calloc(dpamount, sizeof(double*));
    if (dps == NULL) { ErrorOccured(); }
    for (i = 0; i < dpamount ; i++){
        dps[i] = (double*)calloc(dplen, sizeof(double));
        if (dps[i] == NULL) { ErrorOccured(); }
    }
    for (i = 0 ; i < dpamount ; i++){
        for (j = 0; j < dplen ; j++){
            fscanf(input, "%lf%c",&point, &psik);
            dps[i][j] = point;
        }
    }
    fclose(input);
}

double euclideanDistance(double* vectorA, double* vectorB){
    int i;
    double sum_squares = 0;
    for (i = 0 ; i < dplen ; i++){
        sum_squares += pow((vectorA[i]-vectorB[i]),2);
    }
    return sqrt(sum_squares);
}

double** WeightedAdjMat() {
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

double** DiagDegMat(double** weights) {
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

double** main(int argc, char* argv[]) {
    if (argc != 3) {
        InvalidInput();
    }
    char* algo = argv[1];
    char* file = argv[2];
    if ("algo" == "jacobi") {
        
    } else { 
        DataPointsData(file);
        if ("algo" == "wam") {
            WeightedAdjMat();
        } else if ("algo" == "ddg") {
            DiagDegMat(WeightedAdjMat()); /* not good */
        } else if ("algo" == "lnorm") {
            NormalGraphLap();
        } else {
            InvalidInput();
    }

    } return NULL;
}