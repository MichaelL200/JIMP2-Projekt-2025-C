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
#include "output.h"

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
    
    // Inicjalizacja obiektu struktury do metody Lanczosa i obliczeń wartości wektorów własnych macierzy L
    LanczosEigenV lev;
    lanczos_init(&lev, input.v_count, input.v_count);
    // Losowanie dowolnego wektora v₁
    lanczos_v1_init(&lev);
    // Inicjalizacyjny krok iteracyjny metody Lanczosa
    lanczos_initial_step(&lev, L);
    // Inicjalizacja generatora liczb losowych
    srand((unsigned int)time(NULL));
    // Iteracje metody Lanczosa dla j = 2, ..., m
    lanczos(&lev, L);
    // Algorytm QR - obliczanie wartości i wektorów własnych T
    qr_algorithm(&lev);
    // Obliczenie przybliżonych wektorów własnych L
    compute_approximate_eigenvectors(&lev);
    
    // Podział grafu na części
    Result res = clusterization(lev.X, lev.n, config.parts, lev.m, config.margin, A);
    lanczos_free(&lev);
    // Wypisanie wyników klasteryzacji
    print_result(&res);

    // Wypisanie podzielnego grafu do pliku wyjściowego
    write_output(config.output_file, &res, &input, A, input.v_count, config.format);

    // Zwolnienie pamięci
    free_input(&input);
    free(A);
    free(L);

    return EXIT_SUCCESS;
}
