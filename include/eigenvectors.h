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

#endif // EIGENVECTORS_H
