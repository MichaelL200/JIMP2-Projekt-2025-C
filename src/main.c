#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "input.h"
#include "mat_vec.h"

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

    // Wypisanie informacji wczytanych z pliku wejściowego
    print_input(&input);

    // Wczytanie grafu do macierzy sąsiedztwa A
    int *A = get_adjacency_matrix(&input);

    // Obliczenie macierzy stopni grafu D
    int *D = calc_degree_mat(A, input.v_count);

    // Obliczenie macierzy Laplace'a grafu L
    int *L = calc_laplacian(A, D, input.v_count);

    // Zwolnienie pamięci
    free(input.vertices_ptrs);
    free(input.vertices_groups);
    free(A);
    free(D);
    free(L);

    return EXIT_SUCCESS;
}
