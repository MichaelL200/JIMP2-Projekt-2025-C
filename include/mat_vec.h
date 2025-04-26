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

// Obliczenie macierzy stopni grafu D
int* calc_degree_mat(int *A, int n);

// Obliczenie macierzy Laplace'a grafu L
int* calc_laplacian(int* A, int* D, int n);

// Funkcja obliczająca iloczyn skalarny dwóch wektorów o długości n
double dot_product(double *v1, double *v2, int n);

// Funkcja obliczająca normę euklidesową wektora o długości n
double norm(double *v, int n);

// Funkcja mnożąca macierz M (n x n, przechowywana w porządku wierszowym) przez wektor v. Wynik zapisywany w tablicy result
void mat_vec_multiply_d(double* M, double* v, double* result, int n);
void mat_vec_multiply_i(int* M, double* v, double* result, int n);

// Funkcja ortogonalizująca wektor w stosunku do poprzednich wektorów
void orthogonalize(double *v, double *V, int n, int j);

// Funkcja normalizująca wektor
int normalize(double *v, int n);

#endif // MAT_VEC_H
