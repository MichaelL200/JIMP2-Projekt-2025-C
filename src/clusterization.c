#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <limits.h>
#include <omp.h>

#include "clusterization.h"
#include "eigenvectors.h"
#include "output.h"
#include "utils.h"

// Klasteryzacja k-means
int* clusterization(float* eigenvectors, int n, int k, int dim)
{
    // Alokacja pamięci dla tablicy klastrów
    int* clusters = malloc(n * sizeof(int));
    check_alloc(clusters);

    // Alokacja pamięci dla punktów w przestrzeni wektorów własnych
    float** points = malloc(n * sizeof(float*));
    check_alloc(points);
    for (int i = 0; i < n; i++)
    {
        points[i] = malloc(dim * sizeof(float));
        check_alloc(points[i]);
        for (int d = 0; d < dim; d++)
        {
            points[i][d] = eigenvectors[i * dim + d]; // Użycie wektorów własnych bezpośrednio
        }
    }

    // Inicjalizacja centroidów jako pierwszych k punktów
    float** centroids = malloc(k * sizeof(float*));
    check_alloc(centroids);
    for (int i = 0; i < k; i++)
    {
        centroids[i] = malloc(dim * sizeof(float));
        check_alloc(centroids[i]);
        for (int d = 0; d < dim; d++)
        {
            centroids[i][d] = points[i][d];
        }
    }

    // Alokacja pamięci dla liczników i przypisań
    int* counts = malloc(k * sizeof(int));
    check_alloc(counts);
    int* assignments = malloc(n * sizeof(int));
    check_alloc(assignments);
    for (int i = 0; i < n; i++)
    {
        assignments[i] = -1;
    }

    int changed;
    do
    {
        changed = 0;

        // Przypisanie punktów do najbliższego centroidu
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++)
        {
            float min_dist = DBL_MAX;
            int cluster = -1;
            for (int j = 0; j < k; j++)
            {
                float dist = 0.0;
                for (int d = 0; d < dim; d++)
                {
                    float diff = points[i][d] - centroids[j][d];
                    dist += diff * diff;
                }
                if (dist < min_dist)
                {
                    min_dist = dist;
                    cluster = j;
                }
            }
            if (assignments[i] != cluster)
            {
                assignments[i] = cluster;
                changed = 1;
            }
        }

        // Aktualizacja centroidów
        for (int j = 0; j < k; j++)
        {
            for (int d = 0; d < dim; d++)
            {
                centroids[j][d] = 0.0;
            }
            counts[j] = 0;
        }
        #pragma omp parallel for reduction(+:counts[:k])
        for (int i = 0; i < n; i++)
        {
            int c = assignments[i];
            counts[c] += 1;
            for (int d = 0; d < dim; d++)
            {
                centroids[c][d] += points[i][d];
            }
        }
        for (int j = 0; j < k; j++)
        {
            if (counts[j] > 0)
            {
                for (int d = 0; d < dim; d++)
                {
                    centroids[j][d] /= counts[j];
                }
            }
        }
    } while (changed);

    // Przypisanie wyników do tablicy klastrów
    for (int i = 0; i < n; i++)
    {
        clusters[i] = assignments[i];
    }

    // Zwolnienie pamięci
    for (int i = 0; i < n; i++)
    {
        free(points[i]);
    }
    free(points);
    for (int i = 0; i < k; i++)
    {
        free(centroids[i]);
    }
    free(centroids);
    free(counts);
    free(assignments);

    return clusters;
}

// Wypisywanie klastrów (części)
void print_clusters(int* clusters, int n, int k)
{
    printf("\tKlastery (części):\n");
    for (int i = 0; i < k; i++)
    {
        int count = 0;
        printf("\t\tKlaster %d: ", i);
        for (int j = 0; j < n; j++)
        {
            if (clusters[j] == i)
            {
                printf("%d ", j);
                count++;
            }
        }
        printf("\n\t\t\t(%d)\n", count);
    }
}
