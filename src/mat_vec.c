#include <stdlib.h>
#include <math.h>

#include "mat_vec.h"
#include "eigenvectors.h"
#include "input.h"

// Wyświetlenie wektora
void printv(int *v, int n, int n_row)
{
    printf("\t[ ");
    for(int i = 0; i < n; i++)
    {
        if( !( i % n_row) )
        {
            printf("\n\t");
        }
        printf("%d ", v[i]);
    }
    printf("]\n");
}

// Dodanie krawędzi do macierzy sąsiedztwa
void add_edge(int *v, int a, int b, int n)
{
    // [a, b] = [a * n + b]
    v[a * n + b] = 1;
    v[b * n + a] = 1;
    //symetria względem diagonali
}

// Otrzymanie wartości [a, b] wektora
int getv(int *v, int a, int b, int n)
{
    return v[a * n + b];
}

// Wczytanie grafu do macierzy sąsiedztwa A
int* get_adjacency_matrix(Input *i)
{
    int *A = calloc(i->v_count * i->v_count, sizeof(int));
    int p = 0;
    int v = 0;
    if (i->v_count < 30)
    {
        printf("\n\tPołączenia dodane do macierzy sąsiedztwa:");
    }
    for(int it = 0; it < (int)i->g_count; it++)
    {
        // Zaktualizuj wierzchołek
        if(p < (int)i->p_count && it == i->vertices_ptrs[p])
        {
            v = i->vertices_groups[it];
            if (i->v_count < 30)
            {
                printf("\n\t%d - %d: ", it, v);
            }
            p++;
        }
        else
        {
            if (i->v_count < 30)
            {
                printf("\t%d", i->vertices_groups[it]);
            }
            add_edge(A, v, i->vertices_groups[it], i->v_count);
        }
    }
    if (i->v_count < 30)
    {
        printf("\n");
    }

    // Test
    //printf("\]n\t%d\n", getv(A, 93, 90, v_count));

    // Wyświetlenie macierzy sąsiedztwa
    if(i->v_count < 30)
    {
        printf("\n");
        printv(A, i->v_count * i->v_count, i->v_count);
    }

    return A;
}

// Obliczenie macierzy stopni grafu D
int* calc_degree_mat(int *A, int n)
{
    int *D = calloc(n, sizeof(int));
    if(D == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci!\n");
        free(A);
    }

    // dla każdego wierzchołka
    for(int i = 0; i < n; i++)
    {
        // sprawdź liczbę krawędzi w macierzy sąsiedztwa
        for(int j = 0; j < n; j++)
        {
            if(A[i * n + j] != 0)
            {
                // zwiększanie stopnia wierzchołka
                D[i]++;
                
                // zwiększanie stopnia dodatkowo dla pętli
                if(i == j)
                {
                    D[i]++;
                }
            }
        }
    }

    // Wypisanie macierzy (wektora) stopni grafu D (tylko diagonale)
    
    if(n < 30)
    {
        printf("\n");
        printf("\tStopnie wierzchołków:\n");
        printv(D, n, 30);
        printf("\n");
    }

    return D;
}

// Obliczenie macierzy Laplace'a grafu L
int* calc_laplacian(int* A, int* D, int n)
{
    // = 0
    int* L = calloc(n * n, sizeof(int));
    if(L == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci!\n");
        free(A);
        free(D);
    }

    // + D
    for(int i = 0; i < n; i++)
    {
        L[i * n + i] = D[i];
    }

    // - A
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            L[i * n + j] -= A[i * n + j];
        }
    }

    // Wyświetlenie macierzy Laplace'a
    if(n < 30)
    {
        printf("\n");
        printv(L, n * n, n);
    }

    return L;
}

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
