#ifndef EIGENVECTORS_H
#define EIGENVECTORS_H

#include "mat_vec.h"

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
