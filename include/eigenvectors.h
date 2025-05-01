#ifndef EIGENVECTORS_H
#define EIGENVECTORS_H

#include "mat_vec.h"

// Struktura do przechowywania zmiennych potrzebnych do metody Lanczosa oraz do obliczenia przybliżeń wektorów własnych macierzy Laplace'a grafu
typedef struct
{
    // rozmiar macierzy Laplace'a
    int n;
    // liczba iteracji
    int m;
    // baza ortonormalna V z wektorami v (baza przestrzeni Kryłowa)
    double* V;
    // macierz trójdiagonalna T
    double* alpha;
    double* beta;
    // wartości własne macierzy T
    double* theta;
    // wektory własne macierzy T (macierz)
    double* Y;
    // przybliżone wektory własne macierzy L (macierz)
    double* X;
} LanczosEigenV;
// Nazwy takie same jak na Wikipedii - Lanczos algorithm
// https://en.wikipedia.org/wiki/Lanczos_algorithm

// Mnożenie macierzy CSR przez wektor
void csr_matvec(const CSRMatrix_i* A, const double* x, double* y, int n);

// Metoda Lanczosa
void lanczos(const CSRMatrix_i* A, LanczosEigenV* le, int n, int m);

// Wypisywanie wyniku metody Lanczosa
void print_lev(const LanczosEigenV* l);

// Zwalnianie pamięci dla struktury LanczosEigenV
void free_lev(LanczosEigenV* l);

/*
// Tolerancje
#define NORM_TOL 1e-3
#define LANCZOS_TOL 1e-3
#define MAX_ITER 5
#define TOLERANCE 1e-3

// Inicjalizacja wartości obiektu struktury
void lanczos_init(LanczosEigenV *l, int n, int m);

// Losowanie dowolnego wektora v₁
void lanczos_v1_init(LanczosEigenV *l);

// Pierwszy inicjalizacyjny krok iteracyjny metody Lanczosa
void lanczos_initial_step(LanczosEigenV *l, int* A);

// Iteracje metody Lanczosa dla j = 2, ..., m
void lanczos(LanczosEigenV *l, int *A);

// Funkcja budująca macierz T z wartości własnych
double* build_T(LanczosEigenV *l);

// Funkcja sprawdzająca zbieżność macierzy T
int check_convergence(double* T, int m, double tol);

// Algorytm QR z rotacjami Givensa do obliczenia wartości własnych macierzy trójdiagonalnej T
void qr_algorithm(LanczosEigenV *l);

// Obliczenie X = V * Y
void compute_approximate_eigenvectors(LanczosEigenV *l);

// Funkcja zwalniająca pamięć zajmowaną przez strukturę LanczosEigenV
void lanczos_free(LanczosEigenV *l);
*/

#endif // EIGENVECTORS_H
