#ifndef CLUSTERIZATION_H
#define CLUSTERIZATION_H

#include "output.h"
#include "eigenvectors.h"

// Klasteryzacja k-means
int* clusterization(float* eigenvectors, int n, int k, int dim);

// Wypisywanie klastrów (części)
void print_clusters(int* clusters, int n, int k);

#endif // CLUSTERIZATION_H
