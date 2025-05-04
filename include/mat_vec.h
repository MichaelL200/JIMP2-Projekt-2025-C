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
void print_csr_matrix(const CSRMatrix_i* csr, int n);

// Wypisywanie macierzy Laplace'a
void print_laplacian_matrix(CSRMatrix_i* L, int n);

// Zwalnianie pamięci dla struktury CSRMatrix
void free_csr_matrix(CSRMatrix_i* csr);

#endif // MAT_VEC_H
