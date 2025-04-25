#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#include "eigenvectors.h"
#include "mat_vec.h"

// Tolerancje
#define NORM_TOL 1e-10
#define LANCZOS_TOL 1e-12
#define MAX_ITER 1000
#define TOLERANCE 1e-10

// Funkcja obliczająca iloczyn skalarny dwóch wektorów o długości n
double dot_product(double *v1, double *v2, int n)
{
    if(n < 1)
    {
        fprintf(stderr, "n jest mniejsze niż 1!");
        exit(EXIT_FAILURE);
    }

    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += v1[i] * v2[i];
    }
    return sum;
}

// Funkcja obliczająca normę euklidesową wektora o długości n
double norm(double *v, int n)
{
    return sqrt(dot_product(v, v, n));
}

// Funkcja mnożąca macierz M (n x n, przechowywana w porządku wierszowym) przez wektor v
// Wynik zapisywany jest w tablicy result
void mat_vec_multiply_d(double* M, double* v, double* result, int n)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            result[i] += M[i * n + j] * v[j];
        }
    }
}
void mat_vec_multiply_i(int* M, double* v, double* result, int n)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            result[i] += (double)M[i * n + j] * v[j];
        }
    }
}

// Funkcja testująca funkcje: dot_product, norm i mat_vec_multiply
void test1()
{
    printf("\tTEST 1:\n");

    // Testowe dane
    double v1[] = {1.0, 2.0, 3.0};
    double v2[] = {4.0, 5.0, 6.0};
    int n = 3;

    // Test dot_product
    double dp = dot_product(v1, v2, n);
    printf("\tIloczyn skalarny v1 i v2: %f\n", dp);

    // Test norm
    double norm_v1 = norm(v1, n);
    printf("\tNorma v1: %f\n", norm_v1);

    // Test mat_vec_multiply
    double M[] =
    {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    };
    double result[3];
    mat_vec_multiply_d(M, v1, result, n);
    printf("\tMnożenie macierzy M przez wektor v1 daje:\n\t[");
    for (int i = 0; i < n; i++)
    {
        printf("%f", result[i]);
        if (i < n - 1) printf(", ");
    }
    printf("]\n\n");
}

// Funkcja ortogonalizująca wektor w stosunku do poprzednich wektorów
void orthogonalize(double *v, double *V, int n, int j)
{
    for(int k = 0; k < j; ++k)
    {
        double *vk = V + k * n;
        double proj = dot_product(v, vk, n);
        for(int i = 0; i < n; ++i)
        {
            v[i] -= proj * vk[i];
        }
    }
}

// Funkcja normalizująca wektor
int normalize(double *v, int n)
{
    double v_norm = norm(v, n);
    if(v_norm < LANCZOS_TOL)
    {
        return -1; // Norma zbyt mała
    }
    for(int i = 0; i < n; ++i)
    {
        v[i] /= v_norm;
    }
    return 0;
}

// Funkcja testująca funkcje orthogonalize i normalize
void test2()
{
    printf("\tTEST 2:\n");

    int n = 2;

    // Przykład 1: orthogonalizacja [1,1] względem v1=[1,0]
    double V1[2] = { 1.0, 0.0 };    // pojedynczy wektor v1 w V
    double v[2]  = { 1.0, 1.0 };
    printf("\tPrzed orthogonalize: v = [%f, %f]\n", v[0], v[1]);
    orthogonalize(v, V1, n, 1);
    printf("\to orthogonalize:    v = [%f, %f]\n", v[0], v[1]);
    int code = normalize(v, n);
    printf("\tnormalize zwraca %d, v = [%f, %f]\n", code, v[0], v[1]);

    // Przykład 2: normalizacja [3,4]
    double u[2] = { 3.0, 4.0 };
    printf("\tPrzed normalize: u = [%f, %f]\n", u[0], u[1]);
    code = normalize(u, n);
    printf("\tnormalize zwraca %d, u = [%f, %f]\n", code, u[0], u[1]);

    // Przykład 3: próba normalizacji wektora zerowego
    double w[2] = { 0.0, 0.0 };
    printf("\tPrzed normalize: w = [%f, %f]\n", w[0], w[1]);
    code = normalize(w, n);
    printf("\tnormalize zwraca %d (błąd)\n\n", code);
}


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

// Algorytm QR do obliczenia wartości własnych macierzy trójdiagonalnej T
void qr_algorithm(LanczosEigenV *l)
{
    double* T = build_T(l);
    int iter = 0;
    int converged = 0;

    // QR rozkład T = QR (rotacja Givensa)
    while (iter < MAX_ITER && !converged)
    { 
        for (int i = 0; i < l->m - 1; ++i)
        {
            double a = T[i * l->m + i];
            double b = T[(i + 1) * l->m + i];

            double r = hypot(a, b);  // sqrt(a^2 + b^2)
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
        }

        converged = check_convergence(T, l->m, TOLERANCE);
        iter++;
    }

    // Przepisanie wartości własnych z przekątnej macierzy T
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

    // Obliczenie wektorów własnych
    

    free(T);
}

// Rozwiązanie układu równań (T - θᵢ I) x = d (algorytm Thomasa)
void solve_tridiagonal(const double* a, const double* b, const double* c, const double* d, double* x, int n)
{
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
            if (j < n - 1) printf(", ");
        }
        printf("]\n");
    }
}

// Funkcja testująca funkcje metody Lanczosa i algotytm QR
void test3()
{
    printf("\tTEST 3:\n");

    LanczosEigenV l;
    int n = 5;
    // Macierz Laplace'a grafu
    int L[25] =
    {
        2, -1, -1,  0,  0,
       -1,  3, -1, -1,  0,
       -1, -1,  3,  0, -1,
        0, -1,  0,  2, -1,
        0,  0, -1, -1,  2
    };
    printf("\tMacierz Laplace'a grafu:\n");
    for (int i = 0; i < n; i++)
    {
        printf("\t\t"); 
        for (int j = 0; j < n; j++)
        {
            printf("%2d ", L[i * n + j]);
        }
        if(i != n - 1)
        {
            printf("\n");
        }
    }

    lanczos_init(&l, n, n);
    lanczos_v1_init(&l);
    double v1_norm = norm(l.V, n);

    // Wypisanie v1
    printf("\n\tv1 = [ ");
    for(int i = 0; i < n; i++)
    {
        printf("%lf ", l.V[i]);
    }
    printf("]");

    // Test normy
    printf("\tNorma v1: %.6g\n", v1_norm);
    if(fabs(v1_norm - 1.0) < NORM_TOL)
    {
        printf("\tTest udany - v1 znormalizowane!\n");
    }
    else
    {
        printf("\tTest nieudany - v1 nie jest znormalizowane\n");
    }

    // Test niezerowych elementów
    int non_zero_elements = 0;
    for(int i = 0; i < n; i++)
    {
        if(fabs(l.V[i]) > NORM_TOL)
        {
            non_zero_elements++;
        }
    }
    if(non_zero_elements > 0)
    {
        printf("\tTest udany - v1 ma niezerowe elementy");
    }
    else
    {
        printf("Test nieudany - v1 ma wszystkie zerowe elementy");
    }

    lanczos_initial_step(&l, L);

    lanczos(&l, L);

    // Wypisywanie macierzy V
    printf("\tWynikowa macierz V (%d x %d):\n", n, l.m);
    for (int i = 0; i < n; ++i)
    {
        printf("\t\t");
        for (int j = 0; j < l.m; ++j)
        {
            printf("% .3f ", l.V[j*n + i]);
        }
        printf("\n");
    }

    // Wypisywanie macierzy T (m x m)
    // Załóżmy, że już masz l.m, alpha[] i beta[]
    printf("\tWynikowa macierz T (%d x %d):\n", l.m, l.m);
    for (int i = 0; i < l.m; ++i)
    {
        printf("\t\t");
        for (int j = 0; j < l.m; ++j)
        {
            double Tij;
            if (i == j)
            {
                // główna przekątna
                Tij = l.alpha[i];
            }
            else if (j == i + 1)
            {
                // górna przekątna
                Tij = l.beta[j];
            }
            else if (j + 1 == i)
            {
                // dolna przekątna
                Tij = l.beta[i];
            }
            else
            {
                // wszystkie inne elementy zero
                Tij = 0.0;
            }
            printf("% .3f ", Tij);
        }
        if(i != l.m - 1)
        {
            printf("\n");
        }
    }

    // Algorytm QR
    qr_algorithm(&l);

    // Obliczenie wektorów własnych T
    compute_eigenvectors(&l);

    // Obliczenie przybliżonych wektorów własnych L
    compute_approximate_eigenvectors(&l);
    
    free(l.V);
}
