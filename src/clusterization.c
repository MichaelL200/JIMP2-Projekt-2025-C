#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <omp.h>

#include "clusterization.h"
#include "eigenvectors.h"
#include "output.h"
#include "utils.h"

// Comparison function for qsort
int compare_floats(const float* a, const float* b) {
    return (*a > *b) - (*a < *b);
}

// Split clusters based on the Fiedler vector
void split_clusters_by_fiedler_vector(float* fiedler_vector, int n, int* clusters) {
    // Calculate the median of the Fiedler vector
    float* sorted = malloc(n * sizeof(float));
    memcpy(sorted, fiedler_vector, n * sizeof(float));
    qsort(sorted, n, sizeof(float), (int (*)(const void*, const void*))compare_floats);

    float median = sorted[n / 2];
    free(sorted); // Ensure sorted array is freed

    // Assign clusters based on the median
    for (int i = 0; i < n; i++) {
        clusters[i] = (fiedler_vector[i] >= median) ? 1 : 0;
    }
}

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

    // Improved normalization for eigenvectors
    for (int i = 0; i < n; i++) {
        float norm = 0.0;
        for (int d = 0; d < dim; d++) {
            norm += points[i][d] * points[i][d];
        }
        norm = sqrtf(norm);
        if (norm > 0) {
            for (int d = 0; d < dim; d++) {
                points[i][d] /= norm;
            }
        }
    }

    if (k == 2) {
        split_clusters_by_fiedler_vector(eigenvectors + n, n, clusters); // Use the second eigenvector
    } else {
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

        // Enhanced debugging for cluster balance
        float margin = 10.0; // Example margin value
        if (!check_cluster_balance(clusters, n, k, margin)) {
            fprintf(stderr, "Cluster balance check failed. Cluster sizes:\n");
            for (int i = 0; i < k; i++) {
                fprintf(stderr, "Cluster %d: %d points\n", i, counts[i]);
            }
        }

        // Ensure memory for centroids is freed (if k > 2)
        for (int i = 0; i < k; i++) {
            free(centroids[i]); // Free each centroid
        }
        free(centroids); // Free the outer array
        free(counts); // Free counts array
        free(assignments); // Free assignments array
    }

    // Ensure memory for points is freed
    for (int i = 0; i < n; i++) {
        free(points[i]); // Free each inner array
    }
    free(points); // Free the outer array

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

// Sprawdzenie równowagi klastrów
int check_cluster_balance(int* clusters, int n, int k, float margin)
{
    int* counts = malloc(k * sizeof(int));
    check_alloc(counts);
    for (int i = 0; i < k; i++)
    {
        counts[i] = 0;
    }
    for (int i = 0; i < n; i++)
    {
        counts[clusters[i]]++;
    }

    float avg = (float)n / k;
    for (int i = 0; i < k; i++)
    {
        float lower_bound = avg * (1 - margin / 100.0);
        float upper_bound = avg * (1 + margin / 100.0);
        if (counts[i] < lower_bound || counts[i] > upper_bound)
        {
            free(counts); // Ensure counts array is freed
            return 0; // Not balanced
        }
    }

    free(counts); // Ensure counts array is freed
    return 1; // Balanced
}