#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "config.h"
#include "input.h"
#include "mat_vec.h"
#include "eigenvectors.h"
#include "clusterization.h"
#include "output.h"
#include "utils.h"

// Próg dla wielowątkowości
#define MULTITHREADING_THRESHOLD 1

// Maksymalna liczba prób podziału spektralnego
#define MAX_ATTEMPTS 10

// Maksymalny rozmiar grafu do wypisywania danych
int max_print_size = 110;

int main(int argc, char *argv[])
{
    // Parsowanie argumentów
    Config config = parse_args(argc, argv);
    // Weryfikacja pliku wejściowego
    validate_input_file(&config, argc, argv);
    // Wyświetlanie konfiguracji
    print_config(&config);

    // Otwarcie pliku wejściowego i inicjalizacja
    Input input;
    open_init_input(&input, config.input_file);
    // Czytanie pliku linia po linii
    read_input(&input);
    // Sprawdzenie, czy dane wsadowe są poprawne dla tego grafu
    check_input_data(config.parts, input.v_count, &config.margin);
    // Wypisanie informacji wczytanych z pliku wejściowego
    print_input(&input);

    // Ustawienie liczby procesorów
    if(input.v_count > 1000)
    {
        long num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        if (num_threads < 1)
        {
            error("sysconf");
        }
        printf("\n\tLiczba dostępnych rdzeni: %ld\n", num_threads);
        // Ustawienie liczby wątków w OpenMP
        omp_set_num_threads(num_threads);
    }
    else
    {
        omp_set_num_threads(1);
    }

    // Obliczanie macierzy Laplace'a grafu L
    CSRMatrix_i *L = get_laplacian_matrix(&input);
    // Wypisanie macierzy Laplace'a
    print_laplacian_matrix(L, input.v_count);

    // Obliczanie par własnych
    float* eigenvectors;
    float* eigenvalues;
    compute_eigenvectors(L, input.v_count, config.parts, &eigenvectors, &eigenvalues);
    print_eigenpairs(eigenvalues, eigenvectors, config.parts, input.v_count);

    // Klasteryzacja k-means
    int* clusters = clusterization(eigenvectors, input.v_count, config.parts, config.parts);
    print_clusters(clusters, input.v_count, config.parts);
    free(eigenvectors);
    free(eigenvalues);

    // Sprawdzenie równowagi klastrów
    if(check_cluster_balance(clusters, input.v_count, config.parts, config.margin))
    {
        printf("\n\tPodział spektralny zakończony sukcesem!\n");
    }
    else
    {
        printf("\n\tPodział spektralny nie powiódł się! Margines został przekroczony.\n");
        free(clusters);
        free_csr_matrix(L);
        free_input(&input);
        return EXIT_FAILURE; // Exit with failure if balance check fails
    }

    // Wypisanie klastrów do pliku
    //write_output(config.output_file, &input, clusters, input.v_count, config.format);

    free(clusters);
    free_csr_matrix(L);
    free_input(&input);

    return EXIT_SUCCESS;
}
