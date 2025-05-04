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

// Funkcja porównawcza dla qsort
int compare_floats(const float* a, const float* b)
{
    return (*a > *b) - (*a < *b);
}

// Podział klastrów na podstawie wektora Fiedlera
void split_clusters_by_fiedler_vector(float* fiedler_vector, int n, int* clusters)
{
    // Obliczanie mediany wektora Fiedlera
    float* sorted = malloc(n * sizeof(float));
    memcpy(sorted, fiedler_vector, n * sizeof(float));
    qsort(sorted, n, sizeof(float), (int (*)(const void*, const void*))compare_floats);

    float median = sorted[n / 2];
    free(sorted);

    // Przypisanie klastrów na podstawie mediany
    for (int i = 0; i < n; i++) {
        clusters[i] = (fiedler_vector[i] >= median) ? 1 : 0;
    }
}

// Klasteryzacja k-means
int* clusterization(float* eigenvectors, int n, int k, int dim, int margin)
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

    // Normalizacja wektorów własnych
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
        split_clusters_by_fiedler_vector(eigenvectors + n, n, clusters); // Użycie drugiego wektora własnego
    } else {
        // Inicjalizacja centroidów jako pierwszych k punktów
        float** centroids = malloc(k * sizeof(float*));
        check_alloc(centroids);
        for (int i = 0; i < k; i++)
        {
            centroids[i] = malloc(dim * sizeof(float));
            check_alloc(centroids[i]);
            int idx = rand() % n; // Losowy wybór punktu
            for (int d = 0; d < dim; d++)
            {
                centroids[i][d] = points[idx][d];
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

        printf("\n");
        // Debugowanie równowagi klastrów
        int max_margin = 50; // Maksymalny margines (np. 50%)
        while (!check_cluster_balance(clusters, n, k, margin, NULL) && margin <= max_margin)
        {
            fprintf(stderr, "\tUWAGA: MARGINES %.2f%% JEST ZBYT MAŁY DLA OPTYMALNEGO PODZIAŁU.\n", margin);
            margin *= 1.5; // Zwiększenie marginesu przez mnożnik 1.5
            if (margin > max_margin) margin = max_margin; // Ograniczenie do maksymalnego marginesu
            fprintf(stderr, "\tZWIĘKSZONO MARGINES DO: %.2f%% I PONOWIONO PRÓBĘ.\n", margin);
        }

        if (margin >= max_margin) {
            fprintf(stderr, "NIE UDAŁO SIĘ ZNALEŹĆ OPTYMALNEGO MARGINESU W GRANICACH %.2f%%.\n", max_margin);
        }

        for (int i = 0; i < k; i++) {
            free(centroids[i]);
        }
        free(centroids);
        free(counts);
        free(assignments);
    }

    for (int i = 0; i < n; i++) {
        free(points[i]);
    }
    free(points);

    return clusters;
}

// Wypisywanie klastrów (części)
void print_clusters(int* clusters, int n, int k)
{
    if(n < max_print_size)
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
    else
    {
        printf("\n\tGraf jest zbyt duży by wypisać klastry.\n");
    }
}

// Sprawdzenie równowagi klastrów
int check_cluster_balance(int* clusters, int n, int k, float margin, Result* result)
{
    if(result != NULL) result->parts = k;
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
            if (result != NULL) {
                float exceeded_margin = fabs((float)counts[i] - avg) / avg * 100.0;
                result->margin_kept = (int)exceeded_margin; // Zapisanie przekroczonego marginesu
                result->res = 'F';
            }
            free(counts);
            return 0; // Niezrównoważone
        }
    }

    if (result != NULL)
    {
        float max_margin = 0.0f;
        for (int i = 0; i < k; i++)
        {
            float cluster_margin = fabs((float)counts[i] - avg) / avg * 100.0;
            if (cluster_margin > max_margin)
            {
                max_margin = cluster_margin;
            }
        }
        result->margin_kept = (int)max_margin; // Zapisanie rzeczywistego marginesu
        result->res = 'S';
    }

    free(counts);
    return 1; // Zrównoważone
}