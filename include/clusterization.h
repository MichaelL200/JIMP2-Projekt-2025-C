#ifndef CLUSTERIZATION_H
#define CLUSTERIZATION_H

#include "output.h"
#include "eigenvectors.h"

// Algorytm klasteryzacji k-means
int* clusterization(LanczosEigenV* l, EigenvalueIndex* eigvals, int n, int k);

// Sprawdzanie równowagi klastrów w zadanym marginesie procentowym
int check_cluster_balance(int* clusters, int n, int k, int margin_percent);

/*
// Funkcja pomocnicza do liczenia odległości euklidesowej
double euclidean_distance(const double *a, const double *b, int n);

// Algorytm klasteryzacji centroidów (k-means) z minimalizacją liczby przecięć i modyfikacją macierzy sąsiedztwa
Result *clusterization(double *X, int v_count, int parts, int dimensions, double margin_percentage, int *A);

// Wypisanie wyniku klasteryzacji
void print_result(Result *r);
*/

#endif // CLUSTERIZATION_H
