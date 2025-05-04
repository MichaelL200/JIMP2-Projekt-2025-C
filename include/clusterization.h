#ifndef CLUSTERIZATION_H
#define CLUSTERIZATION_H

#include "output.h"
#include "eigenvectors.h"

// Klasteryzacja k-means
int* clusterization(float* eigenvectors, int n, int k, int dim, int margin, Result* result);

// Wypisywanie klastrów (części)
void print_clusters(int* clusters, int n, int k);

// Sprawdzenie równowagi klastrów
int check_cluster_balance(int* clusters, int n, int k, float margin, Result* result);

#endif // CLUSTERIZATION_H
