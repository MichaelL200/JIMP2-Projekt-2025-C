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
    float** V;
    // macierz trójdiagonalna T
    float* alpha;
    float* beta;
    // wartości własne macierzy T
    float* theta;
    // wektory własne macierzy T (macierz)
    float* Y;
    // przybliżone wektory własne macierzy L (macierz)
    float* X;
} LanczosEigenV;
// Nazwy takie same jak na Wikipedii - Lanczos algorithm
// https://en.wikipedia.org/wiki/Lanczos_algorithm

/*
// Tablica struktur zawierających wartość własną oraz jej indeks
typedef struct 
{
    double value;
    int index;
} EigenvalueIndex;
*/

// Poprawione deklaracje w eigenvectors.h
#pragma once
extern void ssaupd_(int *ido, char *bmat, int *n, char *which,
                   int *nev, float *tol, float *resid, int *ncv,
                   float *V, int *ldv, int *iparam, int *ipntr,
                   float *workd, float *workl, int *lworkl, int *info);

extern void sseupd_(int *rvec, char *All, int *select, float *D,
                   float *V_out, int *ldv_out, float *sigma, char *bmat,
                   int *n, char *which, int *nev, float *tol,
                   float *resid, int *ncv, float *V_in, int *ldv_in,
                   int *iparam, int *ipntr, float *workd, float *workl,
                   int *lworkl, int *info);

// Mnożenie macierzy CSR przez wektor
void csr_matvec(const CSRMatrix_i* A, const float* x, float* y, int n);

// Obliczanie wektorów własnych
void compute_eigenvectors(const CSRMatrix_i* graph, int n, int p, float** eigenvectors, float** eigenvalues);

// Wypisanie par własnych
void print_eigenpairs(float *eigenvalues, float *eigenvectors, int p, int n);

/*
// Metoda Lanczosa
void lanczos(const CSRMatrix_i* A, LanczosEigenV* le, int n, int m);

// Wypisywanie wyniku metody Lanczosa
void print_lanczos(const LanczosEigenV* l);

// Funkcja do rotacji Givensa – modyfikuje macierze Q i T
void apply_givens_rotation(float* a, float* b, float* c, float* s);

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
*/

#endif // EIGENVECTORS_H
