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

// Tablica struktur zawierających wartość własną oraz jej indeks
typedef struct 
{
    double value;
    int index;
} EigenvalueIndex;

// Mnożenie macierzy CSR przez wektor
void csr_matvec(const CSRMatrix_i* A, const double* x, double* y, int n);

// Metoda Lanczosa
void lanczos(const CSRMatrix_i* A, LanczosEigenV* le, int n, int m);

// Wypisywanie wyniku metody Lanczosa
void print_lanczos(const LanczosEigenV* l);

// Funkcja do rotacji Givensa – modyfikuje macierze Q i T
void apply_givens_rotation(double* a, double* b, double* c, double* s);

// Funkcja do obliczania wartości i wektorów własnych macierzy trójdiagonalnej T
void qr_algorithm(LanczosEigenV* le);

// Wypisywanie wyniku algorytmu QR
void print_qr(const LanczosEigenV* l);

// Funkcja porównująca
int compare_eigenvalues(const void* a, const void* b);

// Sortowanie wartości i wektorów własnych rosnąco
EigenvalueIndex* sort_eigenvalues(LanczosEigenV* l, int p);

// Zwalnianie pamięci dla struktury LanczosEigenV
void free_lev(LanczosEigenV* l);

#endif // EIGENVECTORS_H
