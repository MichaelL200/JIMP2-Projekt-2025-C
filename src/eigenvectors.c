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

// Metoda Lanczosa
void lanczos(const CSRMatrix_i* A, LanczosEigenV* l, int n, int m)
{
    l->n = n;
    l->m = m;

    // Alokacja pamięci dla V
    l->V = (float**)malloc(m * sizeof(float*));
    check_alloc(l->V);
    for (int i = 0; i < m; ++i)
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
    unsigned int base_seed = (unsigned int)time(NULL) ^ (getpid() << 16);
    #pragma omp parallel
    {
        unsigned int seed = base_seed + (unsigned int)omp_get_thread_num();
        #pragma omp for
        for (int i = 0; i < n; ++i)
        {
            v0[i] = rand_r(&seed) / (float)RAND_MAX;
        }
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
        if (j > 0) {
            float beta_prev = l->beta[j-1];
            #pragma omp parallel for
            for (int i = 0; i < n; ++i)
            {
                w[i] -= alpha * v[i] + beta_prev * v_prev[i];
            }
        }
        else
        {
            #pragma omp parallel for
            for (int i = 0; i < n; ++i)
            {
                w[i] -= alpha * v[i];
            }
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
    float *alpha = l->alpha;
    float *beta = l->beta;

    // Inicjalizacja Y jako macierzy jednostkowej
    l->Y = (float*)malloc(m * m * sizeof(float));
    memset(l->Y, 0, m * m * sizeof(float));
    for (int i = 0; i < m; ++i) l->Y[i*m + i] = 1.0;

    const int max_iter = 100;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Sprawdzenie zbieżności (pomiń dla uproszczenia)
        
        // Rotacje Givensa dla trójdiagonalnej macierzy
        for (int i = 0; i < m-1; ++i) {
            float a = alpha[i];
            float b = beta[i];
            float c, s;
            apply_givens_rotation(&a, &b, &c, &s);
            
            // Aktualizacja alpha i beta
            alpha[i] = c * a - s * b;
            beta[i] = s * a + c * b;
            if (i < m-2) beta[i+1] *= c;

            // Aktualizacja Y
            for (int j = 0; j < m; ++j) {
                float y0 = l->Y[j*m + i];
                float y1 = l->Y[j*m + i+1];
                l->Y[j*m + i]   = c * y0 - s * y1;
                l->Y[j*m + i+1] = s * y0 + c * y1;
            }
        }
    }

    // Oblicz X = V * Y (wektory własne macierzy Laplace'a)
    l->X = (float*)malloc(l->n * l->m * sizeof(float));
    check_alloc(l->X);
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < l->n; ++i) {
        for (int j = 0; j < l->m; ++j) {
            l->X[i * l->m + j] = 0.0;
            for (int k = 0; k < l->m; ++k) {
                l->X[i * l->m + j] += l->V[i][k] * l->Y[k * l->m + j];
            }
        }
    }

    // Przypisz theta z alpha
    l->theta = (float*)malloc(l->m * sizeof(float));
    memcpy(l->theta, alpha, l->m * sizeof(float)); // Kopiuj, nie przypisuj wskaźnika!
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