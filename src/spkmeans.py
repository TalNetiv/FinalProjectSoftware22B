import enum
import sys
import pandas as pd
import numpy as np
import spkmeans
import math

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

def createListForC(points, n, k):
    points = [i.tolist() for i in points]
    c_points = []
    for i in points: #convert all data points to one long list of floats
        for j in i:
            c_points.append(j)
    return c_points 

def reorderEigen(eigens, n):
    mat = [[] for i in range(n+1)]
    j = 0
    while (j < n):
        eigenvals = eigens[0]
        maxi = max(eigenvals)
        maxi_ind = eigenvals.index(maxi)
        for k in range(n+1):
            mat[k][j] = eigens[k][maxi_ind]
        eigenvals[0][maxi_ind] = sys.MIN_SIZE
        j += 1
    return mat

def eigengap(mat, n):
    k = 0
    for i in range(1,n):
        currgap = abs(mat[i-1]-mat[i])
        if currgap > k:
            k = currgap
    return k

def kGreatestCols(eigens, n, k):
    mat = [[0 for m in range(k)] for l in range(n)]
    for i in range(n):
        for j in range(k):
            mat[i][j] = eigens[i][j]
    return mat

def renormalize(U, n, k):
    T = [[] for i in range(n)]
    for i in range(n):
        sum = 0
        for j in range(k):
            sum += (U[i][j])**2
        for j in range(k):
            T[i][j] = (U[i][j])/((sum)**(0.5))
    return T

def initialPoints(data_points,k,N,d):
    np.random.seed(0)
    c = np.random.choice(N) #choose a random first centroid
    initial_indices = [c]
    initial_centroids = np.zeros((1,d)) + data_points[c]
    distance_table = np.full((2,N), math.inf) #initialize a 2 columns, N rows table
    for i in range(k-1):
        distance_table[0] = np.minimum(distance_table[0],np.sum(np.power(data_points - initial_centroids[i],2), axis=1)) #first column of distance table is the distance of each vector from the current centroid
        sumOdds = np.sum(distance_table[0], axis=0)
        np.true_divide(distance_table[0], sumOdds, out=distance_table[1]) #second column of distance table is the probability for each vector to be the next centroid
        c = np.random.choice(N, p=distance_table[1])
        initial_indices.append(c)
        initial_centroids = np.vstack([initial_centroids, data_points[c]])
    return (initial_indices, initial_centroids)


def main(k, goal, filename):
    if goal == goals.JACOBI:
        sym_mat = pd.read_csv()
    else:
        data_points = pd.read_csv(filename, sep = ",", header=None, index_col = 0).sort_index()
        n = len(data_points)
        d = len(data_points.columns)


    if goal == goals.WAM:
        points = createListForC(data_points)
        weightedMat = spkmeans.wam(points, n, d)
        print(weightedMat) #we need to print it in the right way lol

    elif goal == goals.DDG:
        points = createListForC(data_points)
        diagMat = spkmeans.ddg(points, n, d)
        print(diagMat)

    elif goal == goals.LNORM:
        points = createListForC(data_points)
        lapNorm = spkmeans.lnorm(points, n, d)
        print(lapNorm)

    elif goal == goals.JACOBI:
        #how to send the matrix to C?
        pass

    elif goal == goals.SPK:
        try:
            assert (k < n)
        except:
            print("Invalid Input!")
            exit(1)
        mat = createListForC(points)
        lapNorm = spkmeans.lnorm(mat, n, d)
        eigen = spkmeans.jacobi(lapNorm, n)
        orderedEigen = reorderEigen(eigen, n)
        if k == 0:
            k = eigengap(orderedEigen[0], n)
        U = kGreatestCols(orderedEigen[1:], n, k)
        T = renormalize(U, n, k)
        #now we need to create T_points as data points for kmeans somehow
        T_points = T
        (initial_indices, initial_centroids) = initialPoints(T_points, k, n, d)
        c_points = createListForC(T_points)
        c_centroids = createListForC(initial_centroids)
        final_centroids = spkmeans.spk(c_points, c_centroids, n, k, d)
        print(final_centroids, initial_indices)

    else:
        print("Invalid Input!")
        exit(1)
    


checkInput(sys.argv)
main(int(sys.argv[1]), sys.argv[2], sys.argv[3])
    