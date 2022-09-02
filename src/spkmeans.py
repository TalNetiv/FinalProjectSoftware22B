from enum import Enum
import sys
import pandas as pd
import numpy as np
import spkmodule
import math

class goals(Enum):
    SPK = "spk"
    WAM = "wam"
    DDG = "ddg"
    LNORM = "lnorm"
    JACOBI = "jacobi"

def checkInput(input):
    try:
        assert(len(input) == 4)
        assert input[1].isdigit()
    except:
        print("Invalid Input!")
        sys.exit(1)

# converts all data points to one long list of floats
def createListForC(points, numpy=True):
    if numpy:
        points = [i.tolist() for i in points]
    c_points = []
    for i in points:
        for j in i:
            c_points.append(j)
    return c_points

# sort the eigen values returned from jacobi method in decreasing order
def reorderEigen(eigens, n):
    mat = [[0 for m in range(n)] for i in range(n+1)]
    j = 0
    while (j < n):
        eigenvals = eigens[0]
        maxi = max(eigenvals)
        maxi_ind = eigenvals.index(maxi)
        for k in range(n+1):
            mat[k][j] = eigens[k][maxi_ind]
        eigenvals[maxi_ind] = -math.inf
        j += 1
    return mat

# eigengap heuristic function
def eigengap(mat, n):
    max_gap = 0
    k = 0
    for i in range(n//2):
        currgap = abs(mat[i+1]-mat[i])
        if currgap > max_gap:
            max_gap = currgap
            k = i + 1
    return k

# return the k eigen vectors that correspond with the k largest eigen values
def kGreatestCols(eigens, n, k):
    mat = [[0 for m in range(k)] for l in range(n)]
    for i in range(n):
        for j in range(k):
            mat[i][j] = eigens[i][j]
    return mat

def renormalize(U, n, k):
    T = [[0 for m in range(k)] for l in range(n)]
    for i in range(len(U)):
        sum = 0
        for j in range(len(U[0])):
            sum += math.pow((U[i][j]),2)
        if (math.pow(sum, 0.5) == 0):
            continue
        for j in range(len(U[0])):
            T[i][j] = (U[i][j])/(math.pow(sum, 0.5))
    return T

def initialPoints(data_points, k, n, d):
    c = np.random.choice(n-1) #choose a random first centroid
    initial_indices = [c]
    initial_centroids = np.zeros((1,d)) + data_points[c]
    distance_table = np.full((2,n), math.inf) #initialize a 2 columns, N rows table
    for i in range(k-1):
        distance_table[0] = np.minimum(distance_table[0],np.sum(np.power(data_points - initial_centroids[i],2), axis=1)) #first column of distance table is the distance of each vector from the current centroid
        sumOdds = np.sum(distance_table[0], axis=0)
        np.true_divide(distance_table[0], sumOdds, out=distance_table[1]) #second column of distance table is the probability for each vector to be the next centroid
        c = np.random.choice(n, p=distance_table[1])
        initial_indices.append(c)
        initial_centroids = np.vstack([initial_centroids, data_points[c]])
    return (initial_indices, initial_centroids.tolist())

def printMatrix(mat):
    for point in mat:
        print("".join("%.4f," %s for s in point)[:-1])

def printJacobi(mat):
    print("".join("%.4f," %f for f in mat[0]).strip(","))
    for point in mat[1:]:
        print("".join("%.4f," %s for s in point)[:-1])

def printKmeans(initial_indices, final_centroids):
    print("".join("%d," %f for f in initial_indices).strip(","))
    for centroid in final_centroids:
        print("".join("%.4f," %s for s in centroid)[:-1])

def main(k, goal, filename):
    if goal == goals.JACOBI.value:
        sym_mat = pd.read_csv(filename, sep = ",", header=None)
        sym_mat = sym_mat.to_numpy()
        n = len(sym_mat)
    else:
        data_points = pd.read_csv(filename, sep = ",", header=None)
        n = len(data_points)
        d = len(data_points.columns)
        data_points = data_points.to_numpy()

    if goal == goals.WAM.value:
        points = createListForC(data_points)
        weightedMat = spkmodule.wam(points, n, d)
        printMatrix(weightedMat)

    elif goal == goals.DDG.value:
        points = createListForC(data_points)
        diagMat = spkmodule.ddg(points, n, d)
        printMatrix(diagMat)

    elif goal == goals.LNORM.value:
        points = createListForC(data_points)
        lapNorm = spkmodule.lnorm(points, n, d)
        printMatrix(lapNorm)

    elif goal == goals.JACOBI.value:
        points = createListForC(sym_mat)
        jac = spkmodule.jacobi(points, n)
        printJacobi(jac)

    elif goal == goals.SPK.value:
        try:
            assert (k < n and k >= 0)
        except:
            print("Invalid Input!")
            exit(1)
        points = createListForC(data_points)
        lapNorm = spkmodule.lnorm(points, n, d)
        eigens = spkmodule.jacobi(createListForC(lapNorm, False), n)
        orderedEigen = reorderEigen(eigens, n)
        if k == 0:
            k = eigengap(orderedEigen[0], n)
        U = kGreatestCols(orderedEigen[1:], n, k)
        T = renormalize(U, n, k)
        (initial_indices, initial_centroids) = initialPoints(T, k, len(T), len(T[0]))
        if (k == 1):
            printKmeans(initial_indices, initial_centroids)
            exit(0)
        c_points = createListForC(T, False)
        c_centroids = createListForC(initial_centroids, False)
        final_centroids = spkmodule.spk(c_points, c_centroids, len(T), k, len(T[0]))
        printKmeans(initial_indices, final_centroids)
    else:
        print("Invalid Input!")
        exit(1)

checkInput(sys.argv)
np.random.seed(0)
main(int(sys.argv[1]), sys.argv[2], sys.argv[3])
    