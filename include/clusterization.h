#ifndef CLUSTERIZATION_H
#define CLUSTERIZATION_H

// Funkcja pomocnicza do liczenia odległości euklidesowej
double euclidean_distance(const double *a, const double *b, int n);

// Algorytm klasteryzacji centroidów (k-means)
int *clusterization(double *X, int v_count, int parts, int dimensions, double margin_percentage);

#endif // CLUSTERIZATION_H
