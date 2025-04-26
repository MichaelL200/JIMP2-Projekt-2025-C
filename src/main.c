#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "config.h"
#include "input.h"
#include "mat_vec.h"
#include "eigenvectors.h"
#include "test.h"
#include "clusterization.h"

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
    check_input_data(config.parts, input.v_count);
    // Wypisanie informacji wczytanych z pliku wejściowego
    print_input(&input);

    // Wczytanie grafu do macierzy sąsiedztwa A
    int *A = get_adjacency_matrix(&input);
    // Obliczenie macierzy stopni grafu D
    int *D = calc_degree_mat(A, input.v_count);
    // Obliczenie macierzy Laplace'a grafu L
    int *L = calc_laplacian(A, D, input.v_count);
    free(D);

    // Testy poprawności algorytmu
    test1();
    test2();
    test3();
    
    /*
    // Inicjalizacja obiektu struktury do metody Lanczosa i obliczeń wartości wektorów własnych macierzy L
    LanczosEigenV lev;
    lanczos_init(&lev, input.v_count, input.v_count);
    // Losowanie dowolnego wektora v₁
    lanczos_v1_init(&lev);
    lanczos_initial_step(&lev, L);
    srand((unsigned int)time(NULL));
    lanczos(&lev, L);
    qr_algorithm(&lev);
    compute_approximate_eigenvectors(&lev);
    */

    // Zwolnienie pamięci
    free(input.vertices_ptrs);
    free(input.vertices_groups);
    free(A);
    free(L);

    return EXIT_SUCCESS;
}
