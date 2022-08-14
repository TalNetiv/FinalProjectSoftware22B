#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "spkmeans.c"

#ifndef SPKMEANS_H
#define SPKMEANS_H

static void errorOccured();
static double** kMeans(double** points, double** initial_centroids, int n, int k, int d);
static double** weightedAdjMat(double** data_points, int n, int d);
static double** diagDegMat(double** data_points, int n, int d);
static double** normalGraphLap(double** data_points, int n, int d);
static double** jacobian(double** A, int n);

#endif