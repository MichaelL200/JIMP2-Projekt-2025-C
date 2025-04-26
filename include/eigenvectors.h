#ifndef EIGENVECTORS_H
#define EIGENVECTORS_H

// Struktura do przechowywania zmiennych potrzebnych do metody Lanczosa oraz do obliczenia przybliżeń wektorów własnych macierzy Laplace'a grafu
typedef struct
{
    // rozmiar macierzy Laplace'a
    int n;
    // liczba iteracji
    int m;
    // baza ortonormalna V z wektorami v (baza przestrzeni Kryłowa)
    double* V;
    // wektory w (macierz)
    double* W;
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

// Tolerancje
#define NORM_TOL 1e-10
#define LANCZOS_TOL 1e-12
#define MAX_ITER 1000
#define TOLERANCE 1e-10

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

/*
// Rozwiązanie układu równań (T - θᵢ I) x = d (algorytm Thomasa)
void solve_tridiagonal(const double* a, const double* b, const double* c, const double* d, double* x, int n);

// Obliczanie wektorów własnych macierzy T
void compute_eigenvectors(LanczosEigenV *l);
*/

// Obliczenie X = V * Y
void compute_approximate_eigenvectors(LanczosEigenV *l);

#endif // EIGENVECTORS_H
