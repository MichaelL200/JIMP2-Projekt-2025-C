#ifndef MAT_VEC_H
#define MAT_VEC_H

#include "input.h"

// Struktura do przechowywania macierzy w formacie CSR (Compressed Sparse Row)
typedef struct
{
    int nnz;         // Liczba niezerowych elementów
    int* values;  // Wartości niezerowe
    int* col_index;  // Indeksy kolumn
    int* row_ptr;    // Początki wierszy w tablicy values
} CSRMatrix_i;
typedef struct
{
    int nnz;         // Liczba niezerowych elementów
    double* values;  // Wartości niezerowe
    int* col_index;  // Indeksy kolumn
    int* row_ptr;    // Początki wierszy w tablicy values
} CSRMatrix_d;

// Wyświetlenie wektora
void printv(int *v, int n, int n_row);

// Wczytanie grafu do macierzy Laplace'a L
CSRMatrix_i* get_laplacian_matrix(Input* input);

// Wypisywanie macierzy CSR
void print_csr_matrix(CSRMatrix_i* csr, int n);

// Wypisywanie macierzy Laplace'a
void print_laplacian_matrix(CSRMatrix_i* L, int n);

// Zwalnianie pamięci dla struktury CSRMatrix
void free_csr_matrix(CSRMatrix_i* csr);

/*
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
*/

#endif // MAT_VEC_H
