#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "eigenvectors.h"
#include "mat_vec.h"
#include "utils.h"

// Mnożenie macierzy CSR przez wektor
void csr_matvec(const CSRMatrix_i* A, const double* x, double* y, int n)
{
    #pragma omp parallel for
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

// Metoda Lanczosa
void lanczos(const CSRMatrix_i* A, LanczosEigenV* l, int n, int m)
{
    l->n = n;
    l->m = m;

    // Alokacja pamięci
    l->V = (double*)malloc(n * m * sizeof(double));
    check_alloc(l->V);
    l->alpha = (double*)malloc(m * sizeof(double));
    check_alloc(l->alpha);
    l->beta = (double*)malloc((m - 1) * sizeof(double));
    check_alloc(l->beta);
    l->theta = NULL;
    l->Y = NULL;
    l->X = NULL;

    // Inicjalizacja pierwszego wektora v0
    double* v0 = l->V;
    for (int i = 0; i < n; ++i)
    {
        v0[i] = rand() / (double)RAND_MAX;
    }

    // Normalizacja v0
    double norm = 0.0;
    #pragma omp parallel for reduction(+:norm)
    for (int i = 0; i < n; ++i)
    {
        norm += v0[i] * v0[i];
    }
    norm = sqrt(norm);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        v0[i] /= norm;
    }

    // Inicjalizacja zmiennych pomocniczych
    double* v_prev = (double*)calloc(n, sizeof(double));
    check_alloc(v_prev);
    double* w = (double*)malloc(n * sizeof(double));
    check_alloc(w);
    
    for (int j = 0; j < m; ++j)
    {
        double* v = l->V + j * n;

        // w = A * v
        csr_matvec(A, v, w, n);

        // alpha_j = v^T * w
        double alpha = 0.0;
        #pragma omp parallel for reduction(+:alpha)
        for (int i = 0; i < n; ++i)
        {
            alpha += v[i] * w[i];
        }
        l->alpha[j] = alpha;

        // w = w - alpha_j * v - beta_{j-1} * v_{j-1}
        for (int i = 0; i < n; ++i)
        {
            w[i] -= alpha * v[i] + (j > 0 ? l->beta[j - 1] * v_prev[i] : 0.0);
        }

        // beta_j = ||w||
        double beta = 0.0;
        #pragma omp parallel for reduction(+:beta)
        for (int i = 0; i < n; ++i)
        {
            beta += w[i] * w[i];
        }
        beta = sqrt(beta);
        if (j < m - 1)
        {
            l->beta[j] = beta;
        }

        // v_{j+1} = w / beta_j
        if (j < m - 1)
        {
            double* v_next = l->V + (j + 1) * n;
            #pragma omp parallel for
            for (int i = 0; i < n; ++i)
            {
                v_next[i] = w[i] / beta;
            }
        }

        // Aktualizacja v_prev
        memcpy(v_prev, v, n * sizeof(double));
    }

    // Zwolnienie pamięci pomocniczej
    free(v_prev);
    free(w);
}

// Wypisywanie wyniku metody Lanczosa
void print_lanczos(const LanczosEigenV* l)
{
    if(l->n < max_print_size)
    {
        printf("\n\tMacierz trójdiagonalna T:\n");
        printf("\t\tWartości alpha:");
        for (int i = 0; i < l->m; ++i)
        {
            if( !(i % 5) )
            {
                printf("\n\t\t\t");
            }
            printf("%f\t", l->alpha[i]);
        }
        printf("\n\t\tWartości beta:");
        for (int i = 0; i < l->m - 1; ++i)
        {
            if( !(i % 5) )
            {
                printf("\n\t\t\t");
            }
            printf("%f\t", l->beta[i]);
        }

        printf("\n\tMacierz V (wektory ortonormalne):\n");
        for (int i = 0; i < l->n; ++i)
        {
            printf("\t\t");
            for (int j = 0; j < l->m; ++j)
            {
                printf("%10.6f ", l->V[j * l->n + i]);
            }
            printf("\n");
        }
    }
    else
    {
        printf("\n\tGraf jest zbyt duży by wyświetlić wynik metody Lanczosa.\n");
    }
}

// Funkcja do rotacji Givensa – modyfikuje macierze Q i T
void apply_givens_rotation(double* a, double* b, double* c, double* s)
{
    double r = hypot(*a, *b); // pierwiastek z kwadratów a i b
    *c = *a / r;
    *s = -*b / r;
    *a = r;
    *b = 0.0;
}

// Funkcja do obliczania wartości i wektorów własnych macierzy trójdiagonalnej T
void qr_algorithm(LanczosEigenV* l)
{
    int m = l->m;
    int n = l->n;

    // Alokacja pamięci dla macierzy T (m x m)
    double* T = (double*)malloc(m * m * sizeof(double));
    check_alloc(T);
    memset(T, 0, m * m * sizeof(double));
    #pragma omp parallel for
    for (int i = 0; i < m; ++i)
    {
        T[i * m + i] = l->alpha[i];
        if (i < m - 1)
        {
            T[i * m + i + 1] = l->beta[i];
            T[(i + 1) * m + i] = l->beta[i];
        }
    }

    // Alokacja pamięci dla macierzy Y (m x m) i inicjalizacja jako macierz jednostkowa
    l->Y = (double*)malloc(m * m * sizeof(double));
    check_alloc(l->Y);
    memset(l->Y, 0, m * m * sizeof(double));
    for (int i = 0; i < m; ++i)
    {
        l->Y[i * m + i] = 1.0;
    }

    // Iteracje algorytmu QR
    const int max_iter = 1000;
    const double tol = 1e-12; 
    for (int iter = 0; iter < max_iter; ++iter)
    {
        // Sprawdzanie zbieżności
        int converged = 1;
        for (int i = 0; i < m - 1; ++i)
        {
            if (fabs(T[(i + 1) * m + i]) > tol)
            {
                converged = 0;
                break;
            }
        }
        if (converged)
        {
            break;
        }

        // Rotacje Givensa
        for (int i = 0; i < m - 1; ++i)
        {
            double a = T[i * m + i];
            double b = T[(i + 1) * m + i];
            double c, s;
            apply_givens_rotation(&a, &b, &c, &s);

            // Zastosowanie rotacji do T
            for (int j = i; j < m; ++j)
            {
                double tij = T[i * m + j];
                double tji = T[(i + 1) * m + j];
                T[i * m + j] = c * tij - s * tji;
                T[(i + 1) * m + j] = s * tij + c * tji;
            }
            for (int j = 0; j <= i + 1; ++j)
            {
                double tij = T[j * m + i];
                double tji = T[j * m + i + 1];
                T[j * m + i] = c * tij - s * tji;
                T[j * m + i + 1] = s * tij + c * tji;
            }

            // Zastosowanie rotacji do Y (gromadzenie Q)
            for (int j = 0; j < m; ++j)
            {
                double yij = l->Y[j * m + i];
                double yji = l->Y[j * m + i + 1];
                l->Y[j * m + i] = c * yij - s * yji;
                l->Y[j * m + i + 1] = s * yij + c * yji;
            }
        }
    }

    // Po konwergencji, wartości własne znajdują się na diagonali T
    l->theta = (double*)malloc(m * sizeof(double));
    check_alloc(l->theta);
    for (int i = 0; i < m; ++i)
    {
        l->theta[i] = T[i * m + i];
    }

    // Obliczenie przybliżonych wektorów własnych macierzy L: X = V * Y
    l->X = (double*)malloc(n * m * sizeof(double));
    check_alloc(l->X);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            l->X[i * m + j] = 0.0;
            for (int k = 0; k < m; ++k)
            {
                l->X[i * m + j] += l->V[i * m + k] * l->Y[k * m + j];
            }
        }
    }

    free(T);
}

// Wypisywanie wyniku algorytmu QR
void print_qr(const LanczosEigenV* l)
{
    if (l->n < max_print_size)
    {
        // Wypisanie wartości własnych
        printf("\n\tWartości własne:\n");
        for (int i = 0; i < l->m; ++i)
        {
            printf("\t\ttheta_%d = %f\n", i, l->theta[i]);
        }

        // Wypisanie wektorów własnych macierzy T jako wektorów
        printf("\tWektory własne macierzy trójdiagonalnej T:\n");
        for (int j = 0; j < l->m; ++j)
        {
            printf("\t\ty_%d = [", j);
            for (int i = 0; i < l->m; ++i)
            {
                printf("%10.6f", l->Y[i * l->m + j]);
                if (i < l->m - 1) printf(", ");
            }
            printf("]\n");
        }

        // Wypisanie wektorów własnych macierzy Laplace'a L jako wektorów
        printf("\tWektory własne macierzy Laplace'a grafu L:\n");
        for (int j = 0; j < l->m; ++j)
        {
            printf("\t\tx_%d = [", j);
            for (int i = 0; i < l->n; ++i)
            {
                printf("%10.6f", l->X[i * l->m + j]);
                if (i < l->n - 1) printf(", ");
            }
            printf("]\n");
        }

        // Wypisanie 
    }
    else
    {
        printf("\n\tGraf jest zbyt duży, aby wyświetlić wynik algorytmu QR.\n");
    }
}

// Funkcja porównująca
int compare_eigenvalues(const void* a, const void* b)
{
    double diff = ((double*)a)[0] - ((double*)b)[0];
    return (diff > 0) - (diff < 0);
}

// Sortowanie wartości i wektorów własnych rosnąco
EigenvalueIndex *sort_eigenvalues(LanczosEigenV* l, int p)
{
    EigenvalueIndex* eigvals = malloc(l->m * sizeof(EigenvalueIndex));
    check_alloc(eigvals);
    for (int i = 0; i < l->m; ++i)
    {
        eigvals[i].value = l->theta[i];
        eigvals[i].index = i;
    }

    // Sortowanie wartości własnych
    qsort(eigvals, l->m, sizeof(EigenvalueIndex), compare_eigenvalues);

    // Wyświetlanie posortowanych wartości własnych i odpowiadających im wektorów własne
    int count = 0;
    printf("\n\tPosortowane wartości własne:\n");
    for (int i = 0; i < l->m; ++i)
    {
        printf("\t\ttheta_%d = %.6f", eigvals[i].index, eigvals[i].value);
        if(i != 0 && count < p)
        {
            printf("\t[TO CLUSTERIZATION]");
            count++;
        }    
        printf("\n\t\tx_%d = [", eigvals[i].index);
        for (int j = 0; j < l->n; ++j)
        {
            printf(" %.6f", l->X[j * l->m + eigvals[i].index]);
        }
        printf(" ]");
        printf("\n");
    }

    return eigvals;
}

// Zwalnianie pamięci dla struktury LanczosEigenV
void free_lev(LanczosEigenV* l)
{
    if(l == NULL)
    {
        error("Zmienna LanczosEigenV nie może być NULL");
    }

    if(l->V) free(l->V);
    if(l->alpha) free(l->alpha);
    if(l->beta) free(l->beta);
    if(l->theta) free(l->theta);
    if(l->Y) free(l->Y);
    if(l->X) free(l->X);
}