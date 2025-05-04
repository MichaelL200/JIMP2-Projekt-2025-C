#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <limits.h>

#include "clusterization.h"
#include "eigenvectors.h"
#include "output.h"
#include "utils.h"

/*
// Algorytm klasteryzacji k-means
int* clusterization(LanczosEigenV* l, EigenvalueIndex* eigvals, int n, int k)
{
    int dim = k - 1; // Liczba używanych wektorów własnych (pomijając pierwszy)
    int* clusters = malloc(n * sizeof(int));
    check_alloc(clusters);

    // Alokacja przestrzeni dla reprezentacji punktów w przestrzeni własnej
    float** points = malloc(n * sizeof(float*));
    check_alloc(points);
    for (int i = 0; i < n; i++)
    {
        points[i] = malloc(dim * sizeof(float));
        check_alloc(points[i]);
        for (int d = 0; d < dim; d++)
        {
            int eig_idx = eigvals[d + 1].index; // Pomijanie pierwszego wektora własnego
            points[i][d] = l->X[i * l->m + eig_idx];
        }
    }

    // Inicjalizacja centroidów jako pierwsze k punktów
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

    int* counts = malloc(k * sizeof(int));
    check_alloc(counts);
    int* assignments = malloc(n * sizeof(int));
    check_alloc(assignments);
    for (int i = 0; i < n; i++)
        assignments[i] = -1;

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
                centroids[j][d] = 0.0;
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

    // Przypisanie wyników do wektora klastrów
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

// Wypisanie części dla każdego wierzchołka (klastrów)
void print_clusters(int* clusters, int n)
{
    if(n < 3 * max_print_size)
    {
        printf("\n\tCzęści dla każdego wierzchołka:\n");
        printv(clusters, n, 20);
    }
    else
    {
        printf("\n\tGraf jest zbyt duży, by wyświetlić części dla każego wierzchołka.\n");
    }
}

// Sprawdzanie równowagi klastrów w zadanym marginesie procentowym
int check_cluster_balance(int* clusters, int n, int k, int margin)
{
    int* cluster_sizes = calloc(k, sizeof(int));
    check_alloc(cluster_sizes);

    // Zliczanie liczby wierzchołków w każdym klastrze
    for (int i = 0; i < n; i++)
    {
        int cluster_id = clusters[i];
        if (cluster_id < 0 || cluster_id >= k)
        {
            fprintf(stderr, "Nieprawidłowy identyfikator klastra: %d\n", cluster_id);
            free(cluster_sizes);
            return 0;
        }
        cluster_sizes[cluster_id]++;
    }
    free(clusters);

    // Wypisanie rozmiarów klastrów (części)
    printf("\n\tRozmiary części (klastrów):\n");
    for(int i = 0; i < k; i++)
    {
        printf("\t\tCzęść %d: %d\n", i, cluster_sizes[i]);
    }

    // Sprawdzanie najmniejszego i największego klastra (części)
    int min_cluster = n;
    int max_cluster = 0;
    for(int i = 0; i < k; i++)
    {
        if(cluster_sizes[i] < min_cluster)
        {
            min_cluster = cluster_sizes[i];
        }
        if(cluster_sizes[i] > max_cluster)
        {
            max_cluster = cluster_sizes[i];
        }
    }
    free(cluster_sizes);
    int margin_kept = (int)(100.0 * ((double)(max_cluster - min_cluster) / min_cluster));

    // Sprawdzanie marginesu
    if( margin_kept < margin)
    {
        return 1;
    }
    else
    {
        printf("\tMargines nie został spełniony.\n");
        return 0;
    }
}
*/