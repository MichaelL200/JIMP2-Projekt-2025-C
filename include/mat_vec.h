#ifndef MAT_VEC_H
#define MAT_VEC_H

#include "input.h"

// Wyświetlenie wektora
void printv(int *v, int n, int n_row);

// Dodanie krawędzi do macierzy sąsiedztwa
void add_edge(int *v, int a, int b, int n);

// Otrzymanie wartości [a, b] wektora
int getv(int *v, int a, int b, int n);

// Wczytanie grafu do macierzy sąsiedztwa A
int* get_adjacency_matrix(Input *i);

#endif // MAT_VEC_H
