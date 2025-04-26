#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "test.h"
#include "mat_vec.h"
#include "eigenvectors.h"
#include "clusterization.h"

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

// Funkcja testująca funkcje metody Lanczosa i algotytm QR
void test3()
{
    printf("\tTEST 3:\n");

    LanczosEigenV l;
    int n = 5;
    // Macierz sąsiedztwa grafu
    int A[25] =
    {
        0, 1, 1, 0, 0,
        1, 0, 1, 1, 0,
        1, 1, 0, 0, 1,
        0, 1, 0, 0, 1,
        0, 0, 1, 1, 0
    };
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
    printf("\n\tNorma v1: %.6g\n", v1_norm);
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

    // Algorytm QR - obliczanie wartości i wektorów własnych T
    qr_algorithm(&l);

    /*
    // Obliczenie wektorów własnych T
    compute_eigenvectors(&l);
    */

    // Obliczenie przybliżonych wektorów własnych L
    compute_approximate_eigenvectors(&l);

    // Podział grafu na 2 części
    clusterization(l.X, l.n, 2, l.m, 10, A);

    // Wypisanie nowej macierzy sąsiedztwa A
    printf("\tNowa macierz sąsiedztwa A:\n");
    for (int i = 0; i < n; i++)
    {
        printf("\t\t");
        for (int j = 0; j < n; j++)
        {
            printf("%2d ", A[i * n + j]);
        }
        if(i != n - 1)
        {
            printf("\n");
        }
    }
    
    free(l.V);
}
