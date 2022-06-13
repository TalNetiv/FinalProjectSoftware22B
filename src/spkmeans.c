#include <stdio.h>
#include <assert.h>

int dpamount, dplen;
void InvalidInput();
void ErrorOccured();
void DataPointsData(char* name);
int** WeightedAdjMat(char* file);
int** DiagDegMat(char* file);
int** NormalGraphLap(char* file);

void InvalidInput() {
    print('Invalid Input!');
}

void ErrorOccured() {
    print('An Error Has Occured!');
}

void DataPointsData(char* name) {
    FILE * fread = fopen('name', "r");
    dplen = 0, dpamount = 0;
    if (fread == NULL) {
        errorOccured();
        return 1;
    }
    char c = fgetc(fread);
    while (c !='\n'){
        if (c == ','){
            dplen++;
        }
        c = fgetc(fread);
    }
    dplen++, dpamount++;
    while (c != EOF){
        if (c == '\n'){
            dpamount++;
        }
        c = fgetc(fread);
    }
    fclose(fread);
    
}

int** main(int argc, char* argv[]) {
    if (argc != 3) {
        InvalidInput();
    }
    char* algo = argv[1];
    char* file = argv[2];
    if (algo == 'jacobi') {
        
    } else { 
        DataPointsData(file);
        if (algo == 'wam') {
            WeightedAdjMat(file);
        } else if (algo == 'ddg') {
            DiagDegMat(file);
        } else if (algo == 'lnorm') {
            NormalGraphLap(file);
        } else {
            InvalidInput(file);
            return 1;
    }

    } return NULL;
}