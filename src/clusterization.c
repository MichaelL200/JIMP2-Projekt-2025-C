#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

#include "clusterization.h"

#define MAX_ITERATIONS 1000

// Funkcja pomocnicza do liczenia odległości euklidesowej
double euclidean_distance(const double *a, const double *b, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

// Algorytm klasteryzacji centroidów (k-means)
int *clusterization(double *X, int v_count, int parts, int dimensions, double margin_percentage)
{
    if (parts <= 0 || v_count <= 0 || X == NULL || parts > v_count || margin_percentage < 0.0 || margin_percentage > 100.0)
    {
        fprintf(stderr, "Błąd: Niepoprawne dane wejściowe do klasteryzacji.\n");
        exit(EXIT_FAILURE);
    }

    // Alokacja pamięci
    double *centroids = malloc(parts * dimensions * sizeof(double));
    if (!centroids)
    {
        fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla centroidów.\n");
        exit(EXIT_FAILURE);
    }

    int *labels = malloc(v_count * sizeof(int));
    if (!labels)
    {
        free(centroids);
        fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla etykiet.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < v_count; i++)
    {
        labels[i] = -1; // Initialize labels to -1
    }

    // Inicjalizacja centroidów na podstawie pierwszych `parts` punktów z X
    for (int i = 0; i < parts; i++)
    {
        for (int j = 0; j < dimensions; j++)
        {
            centroids[i * dimensions + j] = X[i * dimensions + j];
        }
    }

    // Algorytm k-means z uwzględnieniem ograniczeń
    int iterations = 0;
    int changed;
    int base_vertices_per_cluster = v_count / parts;
    int margin = (int)(base_vertices_per_cluster * (margin_percentage / 100.0));
    // Górny margines liczby wierzchołków w klastrze
    int max_vertices_per_cluster = base_vertices_per_cluster + margin;
    // Dolny margines liczby wierzchołków w klastrze
    int min_vertices_per_cluster = base_vertices_per_cluster - margin;

    // Ensure valid constraints
    if (min_vertices_per_cluster < 0)
    {
        min_vertices_per_cluster = 0;
    }
    if (max_vertices_per_cluster > v_count)
    {
        max_vertices_per_cluster = v_count;
    }

    do
    {
        changed = 0;

        // Przypisanie wierzchołków do najbliższego centroidu z uwzględnieniem ograniczeń
        int *cluster_sizes = calloc(parts, sizeof(int));
        if (!cluster_sizes)
        {
            free(centroids);
            free(labels);
            fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla rozmiarów klastrów.\n");
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < v_count; i++)
        {
            double min_distance = DBL_MAX;
            int best_cluster = -1;

            for (int j = 0; j < parts; j++)
            {
                if (cluster_sizes[j] < max_vertices_per_cluster)
                {
                    double distance = euclidean_distance(&X[i * dimensions], &centroids[j * dimensions], dimensions);
                    if (distance < min_distance)
                    {
                        min_distance = distance;
                        best_cluster = j;
                    }
                }
            }

            if (labels[i] != best_cluster)
            {
                if (labels[i] != -1)
                    cluster_sizes[labels[i]]--; // Zmniejsz licznik starego klastra
                labels[i] = best_cluster;
                cluster_sizes[best_cluster]++; // Zwiększ licznik nowego klastra
                changed = 1;
            }
        }

        // Przeliczenie nowych centroidów
        for (int i = 0; i < parts; i++)
        {
            int count = 0;
            for (int j = 0; j < dimensions; j++)
                centroids[i * dimensions + j] = 0.0;

            for (int j = 0; j < v_count; j++)
            {
                if (labels[j] == i)
                {
                    for (int k = 0; k < dimensions; k++)
                        centroids[i * dimensions + k] += X[j * dimensions + k];
                    count++;
                }
            }

            if (count > 0)
            {
                for (int j = 0; j < dimensions; j++)
                    centroids[i * dimensions + j] /= count;
            }
        }

        free(cluster_sizes);
        iterations++;
    } while (changed && iterations < MAX_ITERATIONS);

    free(centroids);

    // Sprawdzenie, czy przekroczono maksymalną liczbę iteracji
    if (iterations >= MAX_ITERATIONS)
    {
        fprintf(stderr, "Ostrzeżenie: Przekroczono maksymalną liczbę iteracji!\n");
    }

    // Wypisanie wyników klasteryzacji
    printf("\tKlasteryzacja zakończona po %d iteracjach.\n", iterations);
    printf("\tEtykiety klastrów:\n");
    for (int i = 0; i < v_count; i++)
    {
        printf("\t\tWierzchołek %d: Klaster %d\n", i, labels[i]);
    }

    return labels;
}
