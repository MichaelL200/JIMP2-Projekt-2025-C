#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h> // Dołączono dla memset

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

// Wczytanie grafu do macierzy Laplace'a L
CSRMatrix_i* get_laplacian_matrix(Input* input)
{
    int n = input->v_count;
    int* degrees = calloc(n, sizeof(int));
    int* nnz_per_row = calloc(n, sizeof(int));
    int total_nnz = 0;

    // Oblicz stopnie wierzchołków i liczbę niezerowych elementów w każdym wierszu
    for (int i = 0; i < (int)input->p_count; i++)
    {
        int row_start = input->vertices_ptrs[i];
        int row_end = (i + 1 < (int)input->p_count) ? input->vertices_ptrs[i + 1] : (int)input->g_count;
        int vertex = input->vertices_groups[row_start];

        for (int j = row_start + 1; j < row_end; j++) {
            int neighbor = input->vertices_groups[j];
            degrees[vertex]++;
            degrees[neighbor]++; // Zwiększenie stopnia dla sąsiada
            nnz_per_row[vertex]++;
            nnz_per_row[neighbor]++; // Uwzględnienie krawędzi dwukierunkowej
        }
    }

    // Dodanie elementów diagonalnych do nnz_per_row
    for (int i = 0; i < n; i++) {
        nnz_per_row[i]++; // Każdy wiersz ma jeden element diagonalny
        total_nnz += nnz_per_row[i];
    }

    // Alokacja pamięci dla struktury CSR
    CSRMatrix_i* csr = malloc(sizeof(CSRMatrix_i));
    csr->nnz = total_nnz;
    csr->values = malloc(total_nnz * sizeof(double));
    memset(csr->values, 0, total_nnz * sizeof(double)); // Inicjalizacja wartości na zero
    csr->col_index = malloc(total_nnz * sizeof(int));
    memset(csr->col_index, 0, total_nnz * sizeof(int)); // Inicjalizacja col_index na zero
    csr->row_ptr = malloc((n + 1) * sizeof(int));

    if (!csr->values || !csr->col_index || !csr->row_ptr) {
        fprintf(stderr, "Nieudana alokacja pamięci w get_laplacian_matrix.\n");
        free(csr->values); free(csr->col_index); free(csr->row_ptr); free(csr); // Zwolnij pamięć
        free(degrees); free(nnz_per_row); // Zwolnij tymczasowe tablice
        return NULL;
    }

    // Wypełnij row_ptr
    csr->row_ptr[0] = 0;
    for (int i = 0; i < n; i++)
    {
        csr->row_ptr[i + 1] = csr->row_ptr[i] + nnz_per_row[i];
    }

    // Wypełnij values i col_index
    int* current_position = calloc(n, sizeof(int));
    for (int i = 0; i < (int)input->p_count; i++)
    {
        int row_start = input->vertices_ptrs[i];
        int row_end = (i + 1 < (int)input->p_count) ? input->vertices_ptrs[i + 1] : (int)input->g_count;
        int vertex = input->vertices_groups[row_start];

        for (int j = row_start + 1; j < row_end; j++)
        {
            int neighbor = input->vertices_groups[j];

            // Dodaj krawędź dla wierzchołka
            int idx = csr->row_ptr[vertex] + current_position[vertex]++;
            csr->values[idx] = -1.0;
            csr->col_index[idx] = neighbor;

            // Dodaj krawędź dla sąsiada (dwukierunkowa)
            idx = csr->row_ptr[neighbor] + current_position[neighbor]++;
            csr->values[idx] = -1.0;
            csr->col_index[idx] = vertex;
        }
    }

    // Dodaj elementy diagonalne
    for (int i = 0; i < n; i++)
    {
        int idx = csr->row_ptr[i] + current_position[i]++;
        csr->values[idx] = degrees[i];
        csr->col_index[idx] = i;
    }

    // Posortuj kolumny w każdym wierszu
    for (int i = 0; i < n; i++)
    {
        int start = csr->row_ptr[i];
        int end = csr->row_ptr[i + 1];
        // Proste sortowanie bąbelkowe; można zastąpić bardziej efektywnym algorytmem
        for (int j = start; j < end - 1; j++)
        {
            for (int k = start; k < end - 1; k++)
            {
                if (csr->col_index[k] > csr->col_index[k + 1])
                {
                    // Zamień kolumny
                    int temp_col = csr->col_index[k];
                    csr->col_index[k] = csr->col_index[k + 1];
                    csr->col_index[k + 1] = temp_col;

                    // Zamień wartości
                    double temp_val = csr->values[k];
                    csr->values[k] = csr->values[k + 1];
                    csr->values[k + 1] = temp_val;
                }
            }
        }
    }

    free(degrees);
    free(nnz_per_row);
    free(current_position); // Upewnij się, że cała tymczasowa pamięć została zwolniona
    
    return csr;
}

// Wypisywanie macierzy CSR
void print_csr_matrix(const CSRMatrix_i* csr, int n) {
    for (int i = 0; i < n; i++) {
        printf("\t\t");
        for (int j = 0; j < n; j++) {
            double val = 0.0;
            for (int k = csr->row_ptr[i]; k < csr->row_ptr[i + 1]; k++) {
                if (csr->col_index[k] == j) {
                    val = csr->values[k];
                    break;
                }
            }
            printf("%3.0f ", val);
        }
        printf("\n");
    }
}

// Wypisywanie macierzy Laplace'a
void print_laplacian_matrix(CSRMatrix_i* L, int n)
{
    if(n < max_print_size)
    {
        printf("\n\tMacierz Laplace'a grafu L:\n");
        print_csr_matrix(L, n);
    }
}

// Zwalnianie pamięci dla struktury CSRMatrix
void free_csr_matrix(CSRMatrix_i* csr)
{
    if (csr != NULL)
    {
        free(csr->values);
        free(csr->col_index);
        free(csr->row_ptr);
        free(csr);
    }
}
