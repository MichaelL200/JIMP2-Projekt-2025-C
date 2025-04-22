#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#include "eigenvectors.h"
#include "mat_vec.h"

// Tolerancje
#define NORM_TOL 1e-10
#define LANCZOS_TOL 1e-12

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
    printf("]\n");
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
    int n = 2;

    // Przykład 1: orthogonalizacja [1,1] względem v1=[1,0]
    double V1[2] = { 1.0, 0.0 };    // pojedynczy wektor v1 w V
    double v[2]  = { 1.0, 1.0 };
    printf("\n\tPrzed orthogonalize: v = [%f, %f]\n", v[0], v[1]);
    orthogonalize(v, V1, n, 1);
    printf("\to orthogonalize:    v = [%f, %f]\n", v[0], v[1]);
    int code = normalize(v, n);
    printf("\tnormalize zwraca %d, v = [%f, %f]\n\n", code, v[0], v[1]);

    // Przykład 2: normalizacja [3,4]
    double u[2] = { 3.0, 4.0 };
    printf("\tPrzed normalize: u = [%f, %f]\n", u[0], u[1]);
    code = normalize(u, n);
    printf("\tnormalize zwraca %d, u = [%f, %f]\n\n", code, u[0], u[1]);

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

    // Alokacja dla wektorów własnych y – macierz o wymiarach m x m
    l->Y = (double *)calloc(m * m, sizeof(double));
    if (l->Y == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla Y.\n");
        exit(EXIT_FAILURE);
    }

    // Przybliżone wektory własne L (x)
    l->X = NULL;
}

// Losowanie dowolnego wektora v₁
void lanczos_v1_init(LanczosEigenV *l)
{
    // Wypełnianie pierwszego wektora losowymi wartościami
    for (int i = 0; i < l->n; i++)
    {
        l->V[i] = (double)rand() / (double)RAND_MAX;
    }

    // Normalizacja v₁: obliczenie normy i podzielenie każdego elementu przez normę
    double norm_v1 = norm(l->V, l->n);
    if (norm_v1 < NORM_TOL)
    {
        fprintf(stderr, "Błąd: wylosowany wektor ma zbyt małą normę.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < l->n; i++)
    {
        l->V[i] /= norm_v1;
    }
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

    double *w = (double *)malloc(n * sizeof(double));
    if(w == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięcfi dla w.\n");
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

        // Obliczenie wj' = A * vj
        mat_vec_multiply_i(A, vj, w, n);

        // Obliczenie alpha_j = wj' * vj
        alpha[j] = dot_product(w, vj, n);

        // Obliczenie wj = wj' - alpha_j * vj - beta_j * w_{j - 1}
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

// Funkcja testująca funkcje metody Lanczosa
void test3()
{
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
    printf("\n\tMacierz Laplace'a grafu:\n");
    for (int i = 0; i < n; i++)
    {
        printf("\t\t"); 
        for (int j = 0; j < n; j++)
        {
            printf("%2d ", L[i * n + j]);
        }
        printf("\n");
    }
    printf("\n");

    lanczos_init(&l, n, n);
    lanczos_v1_init(&l);
    double v1_norm = norm(l.V, n);

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
        printf("\tTest udany - v1 ma niezerowe elementy\n");
    }
    else
    {
        printf("Test nieudany - v1 ma wszystkie zerowe elementy\n");
    }

    // Wypisanie v1
    printf("\n\tv1 = [ ");
    for(int i = 0; i < n; i++)
    {
        printf("%lf ", l.V[i]);
    }
    printf("]\n");

    lanczos_initial_step(&l, L);

    lanczos(&l, L);

    // Wypisywanie macierzy V
    printf("\n\tWynikowa macierz V (%d x %d):\n", n, l.m);
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
    printf("\n\tWynikowa macierz T (%d x %d):\n", l.m, l.m);
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
        printf("\n");
    }

    free(l.V);
}