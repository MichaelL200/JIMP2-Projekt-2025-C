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
        fprintf(stderr, "n is smaller than 1!");
        exit(EXIT_FAILURE);
    }

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
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
void mat_vec_multiply(double* M, double* v, double* result, int n)
{
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = 0; j < n; j++) {
            result[i] += M[i * n + j] * v[j];
        }
    }
}

// Funkcja testująca poprzednie
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
    double M[] = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    };
    double result[3];
    mat_vec_multiply(M, v1, result, n);
    printf("\tMnożenie macierzy M przez wektor v1 daje: [");
    for (int i = 0; i < n; i++) {
        printf("%f", result[i]);
        if (i < n - 1) printf(", ");
    }
    printf("]\n");
}