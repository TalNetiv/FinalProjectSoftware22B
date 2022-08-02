import enum
import sys
import pandas as pd
import numpy as np
import spkmodule

class goals(enum):
    SPK = "spk"
    WAM = "wam"
    DDG = "ddg"
    LNORM = "lnorm"
    JACOBI = "jacobi"

def checkInput(input):
    try:
        assert(len(input) == 4)
        assert isinstance(input[1], int)
    except:
        print("Invalid Input!")
        sys.exit(1)

# def kmeansppInitial():

# def startJaocbi():
#     Reading file for Jacobi

# def eigengap(filename):
#     1.3

def main(k, goal, filename):
    if goal == goals.JACOBI:
        symm_mat = pd.read_csv()
    else:
        data_points = pd.read_csv(filename, sep = ",", header=None, index_col = 0).sort_index()
        N = len(data_points)
        d = len(data_points.columns)
        try:
            assert (k < N)
        except:
            print("Invalid Input!")
            exit(1)
       # (initial_indices, initial_centroids) = kmeansppInitial(data_points, k, N, d)
    #if k==0:
        #k = eigengap(filename)
        


checkInput(sys.argv)
main(int(sys.argv[1]), sys.argv[2], sys.argv[3])
    