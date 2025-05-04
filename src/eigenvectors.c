#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <cblas.h>
#include <unistd.h>

#include "eigenvectors.h"
#include "mat_vec.h"
#include "utils.h"

// Mnożenie macierzy CSR przez wektor
void csr_matvec(const CSRMatrix_i* A, const float* x, float* y, int n)
{
    #pragma omp parallel for schedule(dynamic, 50)
    for (int i = 0; i < n; ++i)
    {
        y[i] = 0.0;
        int row_start = A->row_ptr[i];
        int row_end = A->row_ptr[i + 1];
        for (int j = row_start; j < row_end; ++j)
        {
            y[i] += A->values[j] * x[A->col_index[j]];
        }
    }
}

// Obliczanie wektorów własnych
void compute_eigenvectors(const CSRMatrix_i* graph, int n, int p, float** eigenvectors, float** eigenvalues)
{
    int ido = 0;
    int info = 0;
    char bmat = 'I';
    char which[] = "SM"; // Compute smallest eigenvalues
    int nev = p;
    float tol = 1e-6;
    float* resid = malloc(n * sizeof(float));
    int ncv = (2 * nev < n) ? 2 * nev : n;
    if (ncv > n) ncv = n;

    float* V = malloc(n * ncv * sizeof(float));
    float* workd = malloc(3 * n * sizeof(float));
    float* workl = malloc(ncv * (ncv + 8) * sizeof(float));
    int lworkl = ncv * (ncv + 8);
    int iparam[11] = {0};
    int ipntr[14] = {0};

    iparam[0] = 1; // Exact shifts
    iparam[2] = 1000; // Max iterations
    iparam[6] = 1; // Mode 1: standard eigenvalue problem

    while (1)
    {
        ssaupd_(&ido, &bmat, &n, which, &nev, &tol, resid,
                &ncv, V, &n, iparam, ipntr, workd, workl,
                &lworkl, &info);

        if (ido == 99) break;

        if (ido == -1 || ido == 1)
        {
            csr_matvec(graph, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1], n);
        }
        else
        {
            fprintf(stderr, "Unsupported ARPACK status: %d\n", ido);
            exit(1);
        }
    }

    if (info != 0)
    {
        fprintf(stderr, "ARPACK error: %d\n", info);
        exit(1);
    }

    int* select = malloc(ncv * sizeof(int));
    memset(select, 0, ncv * sizeof(int));
    check_alloc(select);

    int rvec = 1;
    float* D = malloc(nev * sizeof(float));
    check_alloc(D);
    float sigma = 0.0f;

    sseupd_(&rvec, "A", select, D, V, &n, &sigma, &bmat, &n, which,
            &nev, &tol, resid, &ncv, V, &n, iparam, ipntr,
            workd, workl, &lworkl, &info);

    *eigenvectors = malloc(n * nev * sizeof(float));
    check_alloc(*eigenvectors);
    *eigenvalues = malloc(nev * sizeof(float));
    check_alloc(*eigenvalues);

    // Normalize and store eigenvectors
    for (int i = 0; i < nev; i++)
    {
        (*eigenvalues)[i] = D[i];
        float norm = 0.0f;
        for (int j = 0; j < n; j++)
        {
            float value = V[j * nev + i];
            norm += value * value;
            (*eigenvectors)[i * n + j] = value;
        }
        norm = sqrtf(norm);
        for (int j = 0; j < n; j++)
        {
            (*eigenvectors)[i * n + j] /= norm;
        }
    }

    // Sort eigenvalues and eigenvectors in ascending order
    for (int i = 0; i < nev - 1; i++)
    {
        for (int j = i + 1; j < nev; j++)
        {
            if ((*eigenvalues)[i] > (*eigenvalues)[j])
            {
                // Swap eigenvalues
                float temp_val = (*eigenvalues)[i];
                (*eigenvalues)[i] = (*eigenvalues)[j];
                (*eigenvalues)[j] = temp_val;

                // Swap corresponding eigenvectors
                for (int k = 0; k < n; k++)
                {
                    float temp_vec = (*eigenvectors)[i * n + k];
                    (*eigenvectors)[i * n + k] = (*eigenvectors)[j * n + k];
                    (*eigenvectors)[j * n + k] = temp_vec;
                }
            }
        }
    }

    free(select);
    free(V);
    free(workd);
    free(workl);
    free(D);
    free(resid);
}

// Wypisanie par własnych
void print_eigenpairs(float *eigenvalues, float *eigenvectors, int p, int n)
{
    printf("\n");
    for(int i = 0; i < p; i++)
    {
        printf("\tPair %d: eigenvalue = %.6f\n", i, eigenvalues[i]);
        printf("\t  eigenvector[%d] = [", i);
        for(int j = 0; j < n; j++)
        {
            float value = fabs(eigenvectors[i * n + j]) < 1e-6 ? 0.0f : eigenvectors[i * n + j];
            if(j + 1 < n)
            {
                printf("%.6f, ", value);
            }
            else
            {
                printf("%.6f", value);
            }
        }
        printf("]\n\n");
    }
}
