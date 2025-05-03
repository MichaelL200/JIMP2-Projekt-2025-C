#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <limits.h>

#include "clusterization.h"
#include "eigenvectors.h"
#include "output.h"
#include "utils.h"

// Algorytm klasteryzacji k-means
int* clusterization(LanczosEigenV* l, EigenvalueIndex* eigvals, int n, int k)
{
    int dim = k - 1; // Liczba używanych wektorów własnych (pomijając pierwszy)
    int* clusters = malloc(n * sizeof(int));
    check_alloc(clusters);

    // Alokacja przestrzeni dla reprezentacji punktów w przestrzeni własnej
    double** points = malloc(n * sizeof(double*));
    check_alloc(points);
    for (int i = 0; i < n; i++)
    {
        points[i] = malloc(dim * sizeof(double));
        check_alloc(points[i]);
        for (int d = 0; d < dim; d++)
        {
            int eig_idx = eigvals[d + 1].index; // Pomijanie pierwszego wektora własnego
            points[i][d] = l->X[i * l->m + eig_idx];
        }
    }

    // Inicjalizacja centroidów jako pierwsze k punktów
    double** centroids = malloc(k * sizeof(double*));
    check_alloc(centroids);
    for (int i = 0; i < k; i++)
    {
        centroids[i] = malloc(dim * sizeof(double));
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
        for (int i = 0; i < n; i++)
        {
            double min_dist = DBL_MAX;
            int cluster = -1;
            for (int j = 0; j < k; j++)
            {
                double dist = 0.0;
                for (int d = 0; d < dim; d++)
                {
                    double diff = points[i][d] - centroids[j][d];
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
        for (int i = 0; i < n; i++)
        {
            int cluster = assignments[i];
            for (int d = 0; d < dim; d++)
            {
                centroids[cluster][d] += points[i][d];
            }
            counts[cluster]++;
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

/*
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
Result *clusterization(double *X, int v_count, int parts, int dimensions, double margin_percentage, int *A)
{
    // Sprawdzenie poprawności danych wejściowych
    if (parts <= 0 || v_count <= 0 || X == NULL || A == NULL || parts > v_count || margin_percentage < 0.0 || margin_percentage > 100.0)
    {
        fprintf(stderr, "Błąd: Niepoprawne dane wejściowe do klasteryzacji.\n");
        exit(EXIT_FAILURE);
    }

    int *best_labels = NULL;
    int min_intersection_count = INT_MAX;

    for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++)
    {
        // Alokacja pamięci dla centroidów
        double *centroids = malloc(parts * dimensions * sizeof(double));
        if (!centroids)
        {
            fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla centroidów.\n");
            exit(EXIT_FAILURE);
        }

        // Alokacja pamięci dla etykiet
        int *labels = malloc(v_count * sizeof(int));
        if (!labels)
        {
            free(centroids);
            fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla etykiet.\n");
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < v_count; i++)
        {
            labels[i] = -1; // Inicjalizacja etykiet na -1
        }

        // Inicjalizacja centroidów na podstawie analizy głównych składowych (PCA)
        for (int i = 0; i < parts; i++)
        {
            for (int j = 0; j < dimensions; j++)
            {
                centroids[i * dimensions + j] = 0.0;
            }
        }

        // Obliczanie średniej dla każdego wymiaru
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

            // Alokacja pamięci dla rozmiarów klastrów
            int *cluster_sizes = calloc(parts, sizeof(int));
            if (!cluster_sizes)
            {
                free(centroids);
                free(labels);
                fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla rozmiarów klastrów.\n");
                exit(EXIT_FAILURE);
            }

            // Przypisanie wierzchołków do klastrów
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

            // Aktualizacja centroidów
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

        // Obliczenie liczby przecięć
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

    // Alokacja pamięci dla wyniku
    Result *r = malloc(sizeof(Result));
    if (!r)
    {
        fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla wyniku klasteryzacji.\n");
        exit(EXIT_FAILURE);
    }

    r->cut_count = min_intersection_count;
    int max_cluster_size = 0;
    int min_cluster_size = INT_MAX;

    // Obliczenie rozmiarów klastrów
    int *cluster_sizes = calloc(parts, sizeof(int));
    if (!cluster_sizes)
    {
        fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla rozmiarów części.\n");
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
    r->margin_kept = (int)((double)actual_margin / v_count * 100.0);
    if (r->margin_kept <= margin_percentage)
    {
        r->res = 'S';
    }
    else
    {
        free(best_labels);
        free(r);
        return NULL; // Zwróć NULL, jeśli margines nie został zachowany
    }
    r->parts = parts;

    free(best_labels);
    return r; // Zwróć wskaźnik na wynik
}

// Funkcja do wypisania wyniku klasteryzacji
void print_result(Result *r)
{
    printf("\tWynik klasteryzacji:\n");
    printf("\t\tLiczba przecięć: %d\n", r->cut_count);
    printf("\t\tZachowany margines: %d\n", r->margin_kept);
    printf("\t\tWynik: %c\n", r->res);
}
*/
