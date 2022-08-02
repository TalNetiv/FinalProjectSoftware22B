#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "spkmeans.c"



int main(void) /* test matMultiply */
{
    int i;
    int j;
    int cnt = 1;
    double** multipliedMat;
    double matMultAns[9] = {30.0, 36.0, 42.0, 66.0, 81.0, 96.0, 102.0, 126.0, 150.0};

    double ** mat1 = (double**)calloc(3, sizeof(double*));
    if (mat1 == NULL) { errorOccured(); }
    for (i = 0 ; i < 3 ; i ++) {
        mat1[i] = (double*)calloc(3, sizeof(double));
        if (mat1[i] == NULL) { errorOccured(); }
    }
    for (i = 0 ; i < 3 ; i ++) {
        for (j=0; j< 3; j++){
            mat1[i][j] = cnt;
            cnt ++;
        }
    }
    cnt = 1;
    double ** mat2 = (double**)calloc(3, sizeof(double*));
    if (mat2 == NULL) { errorOccured(); }
    for (i = 0 ; i < 3 ; i ++) {
        mat2[i] = (double*)calloc(3, sizeof(double));
        if (mat2[i] == NULL) { errorOccured(); }
    }
    for (i = 0 ; i < 3 ; i ++) {
        for (j=0; j< 3; j++){
            mat2[i][j] = cnt;
            cnt ++;
        }
    }
    multipliedMat = MatMultiply(mat1, mat2,3);
    cnt = 0;
    for (i = 0; i < 3; i++){
        for (j=0; j<3; j++){
            assert(multipliedMat[i][j] == matMultAns[cnt]);
            cnt ++;
        }
    }
    return 0;
}