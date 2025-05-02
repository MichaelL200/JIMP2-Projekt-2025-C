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

// Maksymalna liczba prób podziału spektralnego
#define MAX_ATTEMPTS 100

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

    // Obliczanie macierzy Laplace'a grafu L
    CSRMatrix_i *L = get_laplacian_matrix(&input);
    // Wypisanie macierzy Laplace'a
    print_laplacian_matrix(L, input.v_count);

    // Metoda Lanczosa
    LanczosEigenV lev;
    srand(time(NULL));
    lanczos(L, &lev, input.v_count, input.v_count);
    print_lanczos(&lev);

    // Algorytm QR
    qr_algorithm(&lev);
    print_qr(&lev);
    
    // Sortowanie wartości i wektorów własnych rosnąco
    EigenvalueIndex* eigenindex = sort_eigenvalues(&lev, config.parts);
    
    /*
    // Licznik prób podziału spektralnego
    int attempts = 0;

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
    */

    free(eigenindex);
    free_lev(&lev);
    free_csr_matrix(L);
    free_input(&input);

    return EXIT_SUCCESS;
}
