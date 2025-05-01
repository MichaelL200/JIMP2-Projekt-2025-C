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

// Maksymalny rozmiar grafu do wypisywania danych
int max_print_size = 110;

/*
#define MAX_ATTEMPTS 100 // Maksymalna liczba prób
*/

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

    // Obliczanie macierzy Laplace'a grafu L
    CSRMatrix_i *L = get_laplacian_matrix(&input);
    // Wypisanie macierzy Laplace'a
    if(input.v_count < max_print_size) print_csr_matrix(L, input.v_count);
    
    /*
    int *D = calc_degree_mat(A, input.v_count);
    // Obliczenie macierzy Laplace'a grafu L
    int *L = calc_laplacian(A, D, input.v_count);
    free(D);
    */

    /*
    // Testy poprawności algorytmu
    test1();
    test2();
    test3();
    */
    
    /*
    int attempts = 0; // Licznik liczby prób

    Result *res = NULL;

    do
    {
        attempts++; // Zwiększenie licznika prób

        // Inicjalizacja struktury LanczosEigenV
        LanczosEigenV lev;
        lanczos_init(&lev, input.v_count, input.v_count);
        lanczos_v1_init(&lev);
        lanczos_initial_step(&lev, L);
        srand((unsigned int)time(NULL));
        lanczos(&lev, L);
        qr_algorithm(&lev);
        compute_approximate_eigenvectors(&lev);

        // Próba klasteryzacji
        res = clusterization(lev.X, lev.n, config.parts, lev.m, config.margin, A);

        // Zwolnienie struktury LanczosEigenV
        lanczos_free(&lev);

        // Sprawdzenie, czy przekroczono maksymalną liczbę prób
        if (attempts >= MAX_ATTEMPTS && res == NULL)
        {
            fprintf(stderr, "Ostrzeżenie: Przekroczono maksymalną liczbę prób (%d). Akceptowanie wyniku z niespełnionym marginesem.\n", MAX_ATTEMPTS);
            break;
        }

    } while (res == NULL); // Powtórz, jeśli wynik klasteryzacji to NULL

    // Sprawdzenie, czy wynik klasteryzacji jest poprawny
    if (res != NULL)
    {
        // Wypisanie wyników klasteryzacji
        print_result(res);

        // Zapisanie wyników do pliku
        write_output(config.output_file, res, &input, A, input.v_count, config.format);

        // Zwolnienie pamięci dla wyniku klasteryzacji
        free(res);
    }
    else
    {
        fprintf(stderr, "Błąd: Wynik klasteryzacji jest NULL. Program zakończył działanie z niespełnionym marginesem.\n");
    }

    // Wypisanie liczby prób
    printf("Liczba prób podziału: %d\n", attempts);

    // Zwolnienie pamięci
    free_input(&input);
    free(L);
    free(A); */

    free_csr_matrix(L);
    free_input(&input);

    return EXIT_SUCCESS;
}
