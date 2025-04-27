#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <limits.h>

#include "clusterization.h"
#include "output.h"

#define MAX_ITERATIONS 100
#define MAX_ATTEMPTS 20

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

// Algorytm klasteryzacji centroidów (k-means) z minimalizacją liczby przecięć i modyfikacją macierzy sąsiedztwa
Result clusterization(double *X, int v_count, int parts, int dimensions, double margin_percentage, int *A)
{
    if (parts <= 0 || v_count <= 0 || X == NULL || A == NULL || parts > v_count || margin_percentage < 0.0 || margin_percentage > 100.0)
    {
        fprintf(stderr, "Błąd: Niepoprawne dane wejściowe do klasteryzacji.\n");
        exit(EXIT_FAILURE);
    }

    int *best_labels = NULL;
    int min_intersection_count = INT_MAX;

    for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++)
    {
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

        // Inicjalizacja centroidów na podstawie analizy głównych składowych (PCA)
        for (int i = 0; i < parts; i++)
        {
            for (int j = 0; j < dimensions; j++)
            {
                centroids[i * dimensions + j] = 0.0;
            }
        }

        // Oblicz średnią dla każdego wymiaru
        double *mean = calloc(dimensions, sizeof(double));
        if (!mean)
        {
            free(centroids);
            free(labels);
            fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla średnich.\n");
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < v_count; i++)
        {
            for (int j = 0; j < dimensions; j++)
            {
                mean[j] += X[i * dimensions + j];
            }
        }
        for (int j = 0; j < dimensions; j++)
        {
            mean[j] /= v_count;
        }

        // Ustawienie centroidów wzdłuż głównych osi danych
        for (int i = 0; i < parts; i++)
        {
            for (int j = 0; j < dimensions; j++)
            {
                centroids[i * dimensions + j] = mean[j] + (i - parts / 2.0) * 0.1; // Rozstawienie centroidów wzdłuż osi
            }
        }

        free(mean);

        // Algorytm k-means z uwzględnieniem ograniczeń

        int iterations = 0;
        int changed;
        int base_vertices_per_cluster = v_count / parts;
        int margin = (int)(base_vertices_per_cluster * (margin_percentage / 100.0));
        int max_vertices_per_cluster = base_vertices_per_cluster + margin;
        int min_vertices_per_cluster = base_vertices_per_cluster - margin;

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
                    {
                        cluster_sizes[labels[i]]--;
                    }
                    labels[i] = best_cluster;
                    cluster_sizes[best_cluster]++;
                    changed = 1;
                }
            }

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
                        {
                            centroids[i * dimensions + k] += X[j * dimensions + k];
                        }
                        count++;
                    }
                }

                if (count > 0)
                {
                    for (int j = 0; j < dimensions; j++)
                    {
                        centroids[i * dimensions + j] /= count;
                    }
                }
            }

            free(cluster_sizes);
            iterations++;
        } while (changed && iterations < MAX_ITERATIONS);

        if (iterations >= MAX_ITERATIONS)
        {
            fprintf(stderr, "Ostrzeżenie: Przekroczono maksymalną liczbę iteracji!\n");
        }

        // Modyfikacja macierzy sąsiedztwa A
        int intersection_count = 0;
        for (int i = 0; i < v_count; i++)
        {
            for (int j = 0; j < v_count; j++)
            {
                if (labels[i] != labels[j])
                {
                    if (A[i * v_count + j] != 0)
                    {
                        intersection_count++;
                    }
                }
            }
        }
        intersection_count /= 2;

        // Sprawdzenie, czy liczba przecięć jest mniejsza niż minimalna
        if (intersection_count < min_intersection_count)
        {
            min_intersection_count = intersection_count;
            free(best_labels);
            best_labels = malloc(v_count * sizeof(int));
            if (!best_labels)
            {
                fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla najlepszych etykiet.\n");
                exit(EXIT_FAILURE);
            }
            for (int i = 0; i < v_count; i++)
            {
                best_labels[i] = labels[i];
            }
        }

        free(labels);
        free(centroids);
    }

    printf("\tEtykiety klastrów:\n");
    for (int i = 0; i < v_count; i++)
    {
        printf("\t\tWierzchołek %d: Klaster %d\n", i, best_labels[i]);
    }

    // Alokacja pamięci dla wyniku
    Result r;
    r.cut_count = min_intersection_count;
    int max_cluster_size = 0;
    int min_cluster_size = INT_MAX;

    // Obliczenie rozmiarów klastrów
    int *cluster_sizes = calloc(parts, sizeof(int));
    if (!cluster_sizes)
    {
        fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla rozmiarów klastrów.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < v_count; i++)
    {
        cluster_sizes[best_labels[i]]++;
    }

    for (int i = 0; i < parts; i++)
    {
        if (cluster_sizes[i] > max_cluster_size)
        {
            max_cluster_size = cluster_sizes[i];
        }
        if (cluster_sizes[i] < min_cluster_size)
        {
            min_cluster_size = cluster_sizes[i];
        }
    }

    free(cluster_sizes);

    // Obliczenie zachowanego marginesu
    int actual_margin = max_cluster_size - min_cluster_size;
    r.margin_kept = (int)((double)actual_margin / v_count * 100.0);
    if (r.margin_kept <= margin_percentage)
    {
        r.res = 'S';
    }
    else
    {
        r.res = 'F';
    }
    r.parts = parts;

    free(best_labels);
    return r;
}

// Wypisanie wyniku klasteryzacji
void print_result(Result *r)
{
    printf("\tWynik klasteryzacji:\n");
    printf("\t\tLiczba przecięć: %d\n", r->cut_count);
    printf("\t\tZachowany margines: %d\n", r->margin_kept);
    printf("\t\tWynik: %c\n", r->res);
}
