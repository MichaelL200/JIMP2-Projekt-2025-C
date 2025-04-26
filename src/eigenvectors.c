#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>

#include "eigenvectors.h"
#include "mat_vec.h"

// Inicjalizacja wartości obiektu struktury
void lanczos_init(LanczosEigenV *l, int n, int m)
{
    // Sprawdzenie poprawności argumentów
    if (l == NULL)
    {
        fprintf(stderr, "Błąd: wskaźnik do struktury LanczosEigenV jest NULL.\n");
        exit(EXIT_FAILURE);
    }
    if (n <= 0)
    {
        fprintf(stderr, "Błąd: rozmiar macierzy n = %d musi być dodatni.\n", n);
        exit(EXIT_FAILURE);
    }
    if (m <= 0)
    {
        fprintf(stderr, "\nBłąd: liczba iteracji m = %d musi być dodatnia.\n", m);
        m = n;
        printf("\tUstawiono liczbę iteracji m = n = %d.\n", n);
    }

    l->n = n;
    l->m = m;

    // Alokacja pamięci dla bazy Kryłowa V: macierz o wymiarach n x m
    l->V = (double *)calloc(n * m, sizeof(double));
    if (l->V == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla V.\n");
        exit(EXIT_FAILURE);
    }

    // Alokacja dla W: macierz n × m
    l->W = (double*) calloc(n * m, sizeof(double));
    if (l->W == NULL) {
        fprintf(stderr, "Błąd alokacji pamięci dla W.\n");
        exit(EXIT_FAILURE);
    }

    // Alokacja dla współczynników alpha (długość m)
    l->alpha = (double *)calloc(m, sizeof(double));
    if (l->alpha == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla alpha.\n");
        exit(EXIT_FAILURE);
    }

    // Alokacja dla współczynników beta (długość m)
    l->beta = (double *)calloc(m, sizeof(double));
    if (l->beta == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla beta.\n");
        exit(EXIT_FAILURE);
    }

    // Alokacja dla wartości własnych theta (długość m)
    l->theta = (double *)calloc(m, sizeof(double));
    if (l->theta == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla theta.\n");
        exit(EXIT_FAILURE);
    }

    // Przybliżone wektory własne T (y) i L (x)
    l->Y = NULL;
    l->X = NULL;
}

// Losowanie dowolnego wektora v₁
void lanczos_v1_init(LanczosEigenV *l)
{
    // Zainicjalizowanie generatora liczb losowych
    srand((unsigned int)time(NULL));

    // Wypełnianie pierwszego wektora losowymi wartościami
    for (int i = 0; i < l->n; i++)
    {
        l->V[i] = (double)rand() / (double)RAND_MAX;
    }

    // Normalizacja v₁
    normalize(l->V, l->n);
}

// Pierwszy inicjalizacyjny krok iteracyjny metody Lanczosa
void lanczos_initial_step(LanczosEigenV *l, int* A)
{
    if (l == NULL)
    {
        fprintf(stderr, "Błąd: wskaźnik do struktury LanczosEigenV jest NULL.\n");
        exit(EXIT_FAILURE);
    }
    
    // Pierwszy wektor v1
    double *v1 = l->V;  // elementy V[0, ..., n-1]
    
    // Alokacja pamięci dla tymczasowego wektora w1'
    double *w1_prime = (double *)malloc(l->n * sizeof(double));
    if (w1_prime == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla w1_prime.\n");
        exit(EXIT_FAILURE);
    }
    
    // w1' = A * v1
    mat_vec_multiply_i(A, v1, w1_prime, l->n);
    
    // Obliczenie α1 = w1' * v1
    double alpha1 = dot_product(w1_prime, v1, l->n);
    l->alpha[0] = alpha1;
    
    // Obliczenie w1 = w1' - α1 * v1
    double *w1 = (double *)malloc(l->n * sizeof(double));
    if (w1 == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla w1.\n");
        free(w1_prime);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < l->n; i++)
    {
        w1[i] = w1_prime[i] - alpha1 * v1[i];
    }
    
    // Wypisanie wyników
    printf("\n\tInicjalizacyjny krok iteracyjny metody Lanczosa:\n");
    printf("\talpha1 = %f\n", alpha1);
    printf("\tw1 = [");
    for (int i = 0; i < l->n; i++)
    {
        if( !(i % 5) )
        {
            printf("\n\t");
        }
        printf("%f", w1[i]);
        if (i < l->n - 1) printf(", ");
    }
    printf(" ]\n");

    // Zapisz w1 do W
    for (int i = 0; i < l->n; i++)
    {
        l->W[i] = w1[i];
    }
    
    l->beta[0] = 0;
    l->alpha[0] = alpha1;

    free(w1_prime);
    free(w1);
}

// Iteracje metody Lanczosa dla j = 2, ..., m
void lanczos(LanczosEigenV *l, int* A)
{
    int n = l->n;
    int m = l->m;
    double *V = l->V;
    double *W = l->W;
    double *alpha = l->alpha;
    double *beta  = l->beta;

    // Inicjalizacja wektora tymczasowego dla wj' a potem dla w
    double *w = (double *)malloc(n * sizeof(double));
    if(w == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla w.\n");
        exit(EXIT_FAILURE);
    }

    // Pętla dla j = 2, ..., m (tu: j = 1, ..., m - 1)
    for (int j = 1; j < m; ++j)
    {
        // Obliczenie beta_j
        double *w_prev = W + (j - 1) * n;
        beta[j] = norm(w_prev, n);

        double *vj = V + j * n;
        double *vj_prev = V + (j - 1) * n;

        // sprawdzenie beta_j
        if(beta[j] > LANCZOS_TOL)
        {
            // zmiana wartości v_j
            for(int i = 0; i < n; ++i)
            {
                vj[i] = w_prev[i] / beta[j];
            }
        }
        else
        {
            // wybór dowolnego wektora v_j o normie 1 ortogonalnego do v1, ..., v{j-1}
            for(int i = 0; i < n; ++i)
            {
                vj[i] = (double)rand() / RAND_MAX;
            }
            orthogonalize(vj, V, n, j);
            printf("\tbeta[%d] = %.20f\n", j, beta[j]);
            if(normalize(vj, n) != 0)
            {
                fprintf(stderr, "Nie udało się znormalizować vj!\n");
                free(w);
                exit(EXIT_FAILURE);
            }
        }

        // Obliczenie wj' = A * vj (w to wj')
        mat_vec_multiply_i(A, vj, w, n);

        // Obliczenie alpha_j = wj' * vj (w to wj')
        alpha[j] = dot_product(w, vj, n);

        // Obliczenie wj = wj' - alpha_j * vj - beta_j * w_{j - 1} (nadpisywanie w)
        for(int i = 0; i < n; ++i)
        {
            w[i] = w[i] - alpha[j] * vj[i] - beta[j] * vj_prev[i];
        }

        // Zapisanie wj do W
        for(int i = 0; i < n; ++i)
        {
            W[j * n + i] = w[i];
        }
    }

    free(w);
}

// Funkcja budująca macierz T z wartości własnych
double* build_T(LanczosEigenV *l)
{
    // Sprawdzenie poprawności argumentów
    if (l == NULL)
    {
        fprintf(stderr, "Błąd: wskaźnik do struktury LanczosEigenV jest NULL.\n");
        exit(EXIT_FAILURE);
    }

    // Inicjalizacaja macierzy T
    double *T = (double *)calloc(l->m * l->m, sizeof(double));
    if (T == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla T.\n");
        exit(EXIT_FAILURE);
    }

    // Zapisanie wartości do macierzy T
    for (int i = 0; i < l->m; ++i)
    {
        for (int j = 0; j < l->m; ++j)
        {
            if (i == j)
            {
                T[i * l->m + j] = l->alpha[i];
            }
            else if (j == i + 1)
            {
                T[i * l->m + j] = l->beta[j];
            }
            else if (j + 1 == i)
            {
                T[i * l->m + j] = l->beta[i];
            }
        }
    }

    return T;
}

// Funkcja sprawdzająca zbieżność macierzy T
int check_convergence(double* T, int m, double tol)
{
    for (int i = 1; i < m; ++i)
    {
        if (fabs(T[i * m + (i - 1)]) > tol)
        {
            return 0; // Nie zbiega
        }
    }

    return 1; // Zbiega
}

// Algorytm QR z rotacjami Givensa do obliczenia wartości własnych macierzy trójdiagonalnej T
void qr_algorithm(LanczosEigenV *l)
{
    // Inicjalizacja macierzy Q jako macierzy jednostkowej
    double* Q_total = (double *)calloc(l->m * l->m, sizeof(double));
    if (Q_total == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla Q_total.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < l->m; ++i)
    {
        Q_total[i * l->m + i] = 1.0;
    }    

    double* T = build_T(l);
    int iter = 0;
    int converged = 0;

    // QR rozkład T = QR (rotacja Givensa)
    while (iter < MAX_ITER && !converged)
    { 
        double* Q = (double *)calloc(l->m * l->m, sizeof(double));
        if (Q == NULL)
        {
            fprintf(stderr, "Błąd alokacji pamięci dla Q.\n");
            free(T);
            free(Q_total);
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < l->m; ++i)
        {
            Q[i * l->m + i] = 1.0;
        }

        // Obliczenie rotacji Givensa
        for (int i = 0; i < l->m - 1; ++i)
        {
            double a = T[i * l->m + i];
            double b = T[(i + 1) * l->m + i];

            double r = hypot(a, b);  // sqrt(a^2 + b^2)
            if (r == 0.0) continue;
            double c = a / r;
            double s = -b / r;

            // Zastosowanie rotacji Givensa z lewej strony: G^T * T
            for (int j = 0; j < l->m; ++j)
            {
                double t1 = c * T[i * l->m + j] - s * T[(i + 1) * l->m + j];
                double t2 = s * T[i * l->m + j] + c * T[(i + 1) * l->m + j];
                T[i * l->m + j] = t1;
                T[(i + 1) * l->m + j] = t2;
            }

            // Zastosowanie rotacji Givensa z prawej strony: T * G
            for (int j = 0; j < l->m; ++j)
            {
                double t1 = c * T[j * l->m + i] - s * T[j * l->m + (i + 1)];
                double t2 = s * T[j * l->m + i] + c * T[j * l->m + (i + 1)];
                T[j * l->m + i] = t1;
                T[j * l->m + (i + 1)] = t2;
            }

            // Zbudowanie macierzy Q z aktualnej rotacji Givensa
            for (int j = 0; j < l->m; ++j)
            {
                double q1 = Q[j * l->m + i];
                double q2 = Q[j * l->m + (i + 1)];
                Q[j * l->m + i] = c * q1 - s * q2;
                Q[j * l->m + (i + 1)] = s * q1 + c * q2;
            }
        }

        // Zaktualizowanie macierzy Q_total (Q_total = Q_total * Q)
        double* Q_tmp = calloc(l->m * l->m, sizeof(double));
        for (int i = 0; i < l->m; ++i)
        {
            for (int j = 0; j < l->m; ++j)
            {
                for (int k = 0; k < l->m; ++k)
                {
                    Q_tmp[i * l->m + j] += Q_total[i * l->m + k] * Q[k * l->m + j];
                }
            }
        }
        memcpy(Q_total, Q_tmp, l->m * l->m * sizeof(double));
        free(Q_tmp);
        free(Q);

        converged = check_convergence(T, l->m, TOLERANCE);
        iter++;
    }

    // Zapisanie wartości własnych z przekątnej macierzy T
    for (int i = 0; i < l->m; ++i)
    {
        l->theta[i] = T[i * l->m + i];
    }

    // Wypisanie wartości własnych
    printf("\n\tWartości własne macierzy T:\n");
    for (int i = 0; i < l->m; ++i)
    {
        printf("\t\ttheta[%d] = %f\n", i, l->theta[i]);
    }
    printf("\tLiczba iteracji: %d\n", iter);
    if (converged)
    {
        printf("\tZbieżność osiągnięta po %d iteracjach.\n", iter);
    }
    else
    {
        printf("\tNie osiągnięto zbieżności po %d iteracjach.\n", iter);
    }

    // Inicjalziacja macierzy Y
    l->Y = (double *)calloc(l->m * l->m, sizeof(double));
    if (l->Y == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla Y.\n");
        free(T);
        free(Q_total);
        exit(EXIT_FAILURE);
    }

    // Zapisanie wektorów własnych macierzy T z Q_total
    for (int i = 0; i < l->m * l->m; ++i)
    {
        l->Y[i] = Q_total[i];
    }

    // Wypisanie wektorów własnych
    printf("\tWektory własne macierzy T:\n");
    for (int i = 0; i < l->m; ++i)
    {
        printf("\t\tY[%d] = [", i);
        for (int j = 0; j < l->m; ++j)
        {
            printf("%f", l->Y[i * l->m + j]);
            if (j < l->m - 1) printf(", ");
        }
        printf("]\n");
    }

    free(T);
}

/*
// Rozwiązanie układu równań (T - θᵢ I) x = d (algorytm Thomasa)
void solve_tridiagonal(const double* a, const double* b, const double* c, const double* d, double* x, int n)
{
    if(n < 1)
    {
        fprintf(stderr, "Błąd: rozmiar macierzy n = %d musi być większy niż 0.\n", n);
        exit(EXIT_FAILURE);
    }

    double* c_prime = malloc(n * sizeof(double));
    double* d_prime = malloc(n * sizeof(double));

    if (c_prime == NULL || d_prime == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci.\n");
        exit(EXIT_FAILURE);
    }

    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n - 1; i++)
    {
        double m = b[i] - a[i - 1] * c_prime[i - 1];
        c_prime[i] = c[i] / m;
        d_prime[i] = (d[i] - a[i - 1] * d_prime[i - 1]) / m;
    }
    // dla ostatniego elementu
    double m = b[n - 1] - a[n - 2] * c_prime[n - 2];
    d_prime[n - 1] = (d[n - 1] - a[n - 2] * d_prime[n - 2]) / m;

    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    free(c_prime);
    free(d_prime);
}

// Obliczanie wektorów własnych macierzy T
void compute_eigenvectors(LanczosEigenV *l)
{
    // Sprawdzenie poprawności argumentów
    if (l == NULL)
    {
        fprintf(stderr, "Błąd: wskaźnik do struktury LanczosEigenV jest NULL.\n");
        exit(EXIT_FAILURE);
    }

    int n = l->n;
    int m = l->m;

    // Alokacja dla wektorów własnych y – macierz o wymiarach m x m
    l->Y = (double *)calloc(m * m, sizeof(double));
    if (l->Y == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla Y.\n");
        exit(EXIT_FAILURE);
    }
    // Alokacja dla wektorów własnych x – macierz o wymiarach n x m
    l->X = (double *)calloc(n * m, sizeof(double));
    if (l->X == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla X.\n");
        exit(EXIT_FAILURE);
    }
    
    // Wektory pomocnicze dla algorytmu Thomas'a
    double *a = calloc(m - 1, sizeof(double)); // podprzekątna
    double *b = calloc(m, sizeof(double));     // przekątna
    double *c = calloc(m - 1, sizeof(double)); // nadprzekątna
    double *d = calloc(m, sizeof(double));     // prawa strona
    double *x = calloc(m, sizeof(double));     // rozwiązanie
    if (!a || !b || !c || !d || !x)
    {
        fprintf(stderr, "Błąd: nie udało się zaalokować pamięci dla wektorów pomocniczych.\n");
        exit(EXIT_FAILURE);
    }

    // Wypełnienie wektorów a, b, c na podstawie macierzy T
    for (int i = 0; i < m; ++i)
    {
        b[i] = l->alpha[i];
        if (i < m - 1)
        {
            a[i] = l->beta[i];
            c[i] = l->beta[i];
        }
    }

    for (int i = 0; i < m; ++i)
    {
        // Przygotowanie prawej strony d (wektor b)
        for (int j = 0; j < m; ++j)
        {
            d[j] = (j == 0) ? 1.0 : 0.0; // wektor jednostkowy jako początkowy
        }

        // Modyfikacja przekątnej dla (T - θᵢ I)
        for (int j = 0; j < m; ++j)
        {
            b[j] = l->alpha[j] - l->theta[i];
        }

        // Rozwiązanie układu równań (T - θᵢ I) x = d (algorytm Thomasa)
        solve_tridiagonal(a, b, c, d, x, m);

        // Normalizacja wektora x
        double norm = 0.0;
        for (int j = 0; j < m; ++j)
        {
            norm += x[j] * x[j];
        }
        norm = sqrt(norm);
        for (int j = 0; j < m; ++j)
        {
            l->Y[j * m + i] = x[j] / norm;
        }
    }

    free(a);
    free(b);
    free(c);
    free(d);
    free(x);

    // Wypisanie wektorów własnych T
    printf("\tWektory własne macierzy T:\n");
    for (int i = 0; i < m; ++i)
    {
        printf("\t\tY[%d] = [", i);
        for (int j = 0; j < m; ++j)
        {
            printf("%f", l->Y[j * m + i]);
            if (j < m - 1) printf(", ");
        }
        printf("]\n");
    }
}
*/

// Obliczenie X = V * Y
void compute_approximate_eigenvectors(LanczosEigenV *l)
{
    int n = l->n;
    int m = l->m;
    l->X = calloc(n * m, sizeof(double));
    if (l->X == NULL)
    {
        fprintf(stderr, "Błąd: nie udało się zaalokować pamięci dla macierzy X.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < m; ++k)
            {
                sum += l->V[j * m + k] * l->Y[k * m + i];
            }
            l->X[j * m + i] = sum;
        }
    }

    // Wypisanie przybliżonych wektorów własnych L
    printf("\tPrzybliżone wektory własne L:\n");
    for (int i = 0; i < m; ++i)
    {
        printf("\t\tX[%d] = [", i);
        for (int j = 0; j < n; ++j)
        {
            printf("%f", l->X[j * m + i]);
            if (j < n - 1)
            {
                printf(", ");
            }
        }
        printf("]\n");
    }
}
