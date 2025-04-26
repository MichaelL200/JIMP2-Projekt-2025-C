#ifndef CLUSTERIZATION_H
#define CLUSTERIZATION_H

// Funkcja pomocnicza do liczenia odległości euklidesowej
double euclidean_distance(const double *a, const double *b, int n);

// Algorytm klasteryzacji centroidów (k-means) z minimalizacją liczby przecięć i modyfikacją macierzy sąsiedztwa
int *clusterization(double *X, int v_count, int parts, int dimensions, double margin_percentage, int *A);

#endif // CLUSTERIZATION_H
