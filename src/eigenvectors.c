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
void csr_matvec(const CSRMatrix_i* A, const float* x, float* y, int n)
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

    // Alokacja pamięci dla V
    l->V = (float**)malloc(n * m * sizeof(float));
    check_alloc(l->V);
    for (int i = 0; i < m; ++i) // dla każdego wektora w tablicy V
    {
        l->V[i] = (float*)malloc(n * sizeof(float));
        check_alloc(l->V[i]);
    }
    l->alpha = (float*)malloc(m * sizeof(float));
    check_alloc(l->alpha);
    l->beta = (float*)malloc((m - 1) * sizeof(float));
    check_alloc(l->beta);
    l->theta = NULL;
    l->Y = NULL;
    l->X = NULL;

    // Inicjalizacja pierwszego wektora v0
    float* v0 = l->V[0];
    for (int i = 0; i < n; ++i)
    {
        v0[i] = rand() / (float)RAND_MAX;
    }

    // Normalizacja v0
    float norm = 0.0;
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
    float* v_prev = (float*)calloc(n, sizeof(float));
    check_alloc(v_prev);
    float* w = (float*)malloc(n * sizeof(float));
    check_alloc(w);
    
    for (int j = 0; j < m; ++j)
    {
        float* v = l->V[j];  // Dostęp do j-tego wektora

        // w = A * v
        csr_matvec(A, v, w, n);

        // alpha_j = v^T * w
        float alpha = 0.0;
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
        float beta = 0.0;
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
            float* v_next = l->V[j + 1];  // Dostęp do (j+1)-tego wektora
            #pragma omp parallel for
            for (int i = 0; i < n; ++i)
            {
                v_next[i] = w[i] / beta;
            }
        }

        // Aktualizacja v_prev
        memcpy(v_prev, v, n * sizeof(float));
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
                printf("%10.6f ", l->V[i][j]);
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
void apply_givens_rotation(float* a, float* b, float* c, float* s)
{
    float r = hypot(*a, *b); // pierwiastek z kwadratów a i b
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
    float* T = (float*)malloc(m * m * sizeof(float));
    check_alloc(T);
    memset(T, 0, m * m * sizeof(float));
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
    l->Y = (float*)malloc(m * m * sizeof(float));
    check_alloc(l->Y);
    memset(l->Y, 0, m * m * sizeof(float));
    for (int i = 0; i < m; ++i)
    {
        l->Y[i * m + i] = 1.0;
    }

    // Iteracje algorytmu QR
    const int max_iter = 1000;
    const float tol = 1e-12; 
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
            float a = T[i * m + i];
            float b = T[(i + 1) * m + i];
            float c, s;
            apply_givens_rotation(&a, &b, &c, &s);

            // Zastosowanie rotacji do T
            for (int j = i; j < m; ++j)
            {
                float tij = T[i * m + j];
                float tji = T[(i + 1) * m + j];
                T[i * m + j] = c * tij - s * tji;
                T[(i + 1) * m + j] = s * tij + c * tji;
            }
            for (int j = 0; j <= i + 1; ++j)
            {
                float tij = T[j * m + i];
                float tji = T[j * m + i + 1];
                T[j * m + i] = c * tij - s * tji;
                T[j * m + i + 1] = s * tij + c * tji;
            }

            // Zastosowanie rotacji do Y (gromadzenie Q)
            for (int j = 0; j < m; ++j)
            {
                float yij = l->Y[j * m + i];
                float yji = l->Y[j * m + i + 1];
                l->Y[j * m + i] = c * yij - s * yji;
                l->Y[j * m + i + 1] = s * yij + c * yji;
            }
        }
    }

    // Po konwergencji, wartości własne znajdują się na diagonali T
    l->theta = (float*)malloc(m * sizeof(float));
    check_alloc(l->theta);
    for (int i = 0; i < m; ++i)
    {
        l->theta[i] = T[i * m + i];
    }

    // Obliczenie przybliżonych wektorów własnych macierzy L: X = V * Y
    l->X = (float*)malloc(n * m * sizeof(float));
    check_alloc(l->X);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            l->X[i * m + j] = 0.0;
            for (int k = 0; k < m; ++k)
            {
                l->X[i * m + j] += l->V[i][k] * l->Y[k * m + j];
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

    // Wyświetlanie posortowanych wartości własnych i odpowiadających im wektorów własnych
    int count = 0;
    printf("\n\tPosortowane wartości własne:\n");
    for (int i = 0; i < l->m; ++i)
    {
        printf("\t\ttheta_%d = %.6f", eigvals[i].index, eigvals[i].value);
        if(i != 0 && count < p)
        {
            printf("\t[TO CLUSTERIZATION]"); // Zawiera informacje do klasteryzacji
            count++;
        }    
        printf("\n\t\tx_%d = [", eigvals[i].index);
        for (int j = 0; j < l->n; ++j)
        {
            printf(" %.6f", l->V[eigvals[i].index][j]);  // Pamiętaj, że zmieniłeś strukturę na tablicę tablic
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

    if(l->V)
    {
        for (int i = 0; i < l->m; i++)
        {
            free(l->V[i]);
        }
        free(l->V);
    }
    if(l->alpha) free(l->alpha);
    if(l->beta) free(l->beta);
    if(l->theta) free(l->theta);
    if(l->Y) free(l->Y);
    if(l->X) free(l->X);
}