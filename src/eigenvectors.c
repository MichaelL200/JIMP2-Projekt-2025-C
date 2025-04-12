#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#include "eigenvectors.h"
#include "mat_vec.h"

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
            result[i] += M[i * n + j] * v[j];
        }
    }
}

// Funkcja testująca funkcje: dot_product, norm i mat_vec_multiply
void test_ev()
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

    // Ustawienie pól rozmiaru
    l->n = n;
    l->m = m;

    // Alokacja pamięci dla bazy Q: macierz o wymiarach n x (m+1)
    l->Q = (double *)calloc(n * (m + 1), sizeof(double));
    if (l->Q == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla Q.\n");
        exit(EXIT_FAILURE);
    }

    // Alokacja dla współczynników alpha (długość m)
    l->alpha = (double *)calloc(m, sizeof(double));
    if (l->alpha == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla alpha.\n");
        exit(EXIT_FAILURE);
    }

    // Alokacja dla współczynników beta (długość m+1; beta[0] = 0)
    l->beta = (double *)calloc(m + 1, sizeof(double));
    if (l->beta == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla beta.\n");
        exit(EXIT_FAILURE);
    }

    // Alokacja dla wartości własnych (theta) – długość m
    l->theta = (double *)calloc(m, sizeof(double));
    if (l->theta == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla theta.\n");
        exit(EXIT_FAILURE);
    }

    // Alokacja dla wektorów własnych y – macierz o wymiarach m x m
    l->y = (double *)calloc(m * m, sizeof(double));
    if (l->y == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla y.\n");
        exit(EXIT_FAILURE);
    }

    // Przybliżone wektory własne L (x)
    l->x = NULL;
}

// Losowanie dowolnego wektora v₁
void lanczos_v1_init(LanczosEigenV *l)
{
    // Inicjalizacja pierwszego wektora v₁ ∈ ℂⁿ
    srand((unsigned int)time(NULL));

    // Wypełnianie pierwszego wektora losowymi wartościami
    for (int i = 0; i < l->n; i++)
    {
        l->Q[i] = (double)rand() / (double)RAND_MAX;
    }

    // Normalizacja v₁: obliczenie normy i podzielenie każdego elementu przez normę
    double sum_sq = 0.0;
    for (int i = 0; i < l->n; i++)
    {
        sum_sq += l->Q[i] * l->Q[i];
    }
    double norm_v1 = sqrt(sum_sq);
    if (norm_v1 < 1e-10)
    {
        fprintf(stderr, "Błąd: wylosowany wektor ma zbyt małą normę.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < l->n; i++)
    {
        l->Q[i] /= norm_v1;
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
    double *v1 = l->Q;  // elementy Q[0, ..., n-1]
    
    // Alokacja pamięci dla tymczasowego wektora w1'
    double *w1_prime = (double *)malloc(l->n * sizeof(double));
    if (w1_prime == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla w1_prime.\n");
        exit(EXIT_FAILURE);
    }
    
    // w1' = A * v1
    mat_vec_multiply_i(A, v1, w1_prime, l->n);
    
    // Obliczanie α1
    double alpha1 = dot_product(v1, w1_prime, l->n);
    l->alpha[0] = alpha1;
    
    // Obliczanie w1 = w1' - α1 * v1
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
    
    free(w1_prime);
    free(w1);
}