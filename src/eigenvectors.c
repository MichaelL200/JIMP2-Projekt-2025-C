#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <cblas.h>
#include <unistd.h>

#include "eigenvectors.h"
#include "mat_vec.h"
#include "utils.h"

// Mnożenie macierzy CSR przez wektor
void csr_matvec(const CSRMatrix_i* A, const float* x, float* y, int n)
{
    for (int i = 0; i < n; ++i)
    {
        y[i] = 0.0;
        int row_start = A->row_ptr[i];
        int row_end = A->row_ptr[i + 1];
        for (int j = row_start; j < row_end; ++j)
        {
            y[i] += A->values[j] * x[A->col_index[j]];
        }
    }
}

// Obliczanie wektorów własnych
void compute_eigenvectors(const CSRMatrix_i* graph, int n, int p, float** eigenvectors, float** eigenvalues)
{
    // Sprawdzenie poprawności danych wejściowych
    // n: liczba wierzchołków grafu (rozmiar macierzy)
    // p: liczba par własnych do obliczenia
    if (n <= 0 || p <= 0 || p > n)
    {
        fprintf(stderr, "Nieprawidłowe dane wejściowe: n = %d, p = %d. Upewnij się, że n > 0 i 0 < p <= n.\n", n, p);
        exit(EXIT_FAILURE);
    }

    // Wyświetlenie informacji o rozpoczęciu obliczeń
    printf("Obliczanie %d. par własnych...\n", p);

    // ARPACK: Inicjalizacja zmiennych dla algorytmu iteracyjnego
    int ido = 0; // Status algorytmu (kontrola iteracji)
    int info = 0; // Kod błędu zwracany przez ARPACK
    char bmat = 'I'; // Typ macierzy: 'I' oznacza macierz jednostkową (standardowy problem wartości własnych)
    char which[] = "SM"; // Kryterium wyboru: "SM" oznacza najmniejsze wartości własne (Smallest Magnitude)
    int nev = p; // Liczba wartości własnych do obliczenia
    float tol = 1e-6; // Tolerancja błędu (dokładność obliczeń)

    // Wektor rezydualny (przechowuje aktualne przybliżenie wartości własnych)
    float* resid = malloc(n * sizeof(float));

    // Liczba wektorów Arnoldiego (ncv): wpływa na stabilność i szybkość algorytmu
    int ncv = (4 * nev < n) ? 4 * nev : n;
    if (ncv > n) ncv = n;

    // Alokacja pamięci dla wektorów Arnoldiego i innych struktur roboczych
    float* V = malloc(n * ncv * sizeof(float)); // Macierz wektorów Arnoldiego
    float* workd = malloc(3 * n * sizeof(float)); // Wektory robocze dla ARPACK
    float* workl = malloc(ncv * (ncv + 8) * sizeof(float)); // Tablica robocza dla ARPACK
    int lworkl = ncv * (ncv + 8); // Rozmiar tablicy workl
    int iparam[11] = {0}; // Parametry sterujące ARPACK
    int ipntr[14] = {0}; // Wskaźniki do tablic roboczych

    // Konfiguracja parametrów ARPACK
    iparam[0] = 1; // Dokładne przesunięcia
    iparam[2] = 10000; // Maksymalna liczba iteracji
    iparam[6] = 1; // Tryb 1: standardowy problem wartości własnych

    // Algorytm iteracyjny ARPACK (ssaupd_)
    int retries = 3; // Liczba prób w przypadku niepowodzenia
    while (retries > 0)
    {
        while (1)
        {
            // Wywołanie ARPACK: ssaupd_ wykonuje iteracje w celu obliczenia wartości własnych
            ssaupd_(&ido, &bmat, &n, which, &nev, &tol, resid,
                    &ncv, V, &n, iparam, ipntr, workd, workl,
                    &lworkl, &info);

            // Sprawdzenie błędów
            if (info < 0)
            {
                fprintf(stderr, "Błąd ARPACK ssaupd_: %d. Sprawdź parametry wejściowe.\n", info);
                free(resid); free(V); free(workd); free(workl);
                exit(EXIT_FAILURE);
            }

            // Jeśli ido == 99, oznacza to zakończenie iteracji
            if (ido == 99) break;

            // Obsługa mnożenia macierzy przez wektor (zgodnie z żądaniem ARPACK)
            if (ido == -1 || ido == 1)
            {
                csr_matvec(graph, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1], n);
            }
            else
            {
                fprintf(stderr, "Nieobsługiwany status ARPACK: %d\n", ido);
                free(resid); free(V); free(workd); free(workl);
                exit(1);
            }
        }

        // Sprawdzenie, czy algorytm zakończył się sukcesem
        if (info == 0) break;

        // Jeśli dokładność nie została osiągnięta, zmniejsz tolerancję i spróbuj ponownie
        if (info == 1)
        {
            fprintf(stderr, "Błąd ARPACK: Nieosiągnięta dokładność. Zmniejszanie tolerancji i ponowna próba...\n");
            tol /= 10;
            retries--;
        }
        else
        {
            fprintf(stderr, "Błąd ARPACK: %d\n", info);
            free(resid); free(V); free(workd); free(workl);
            exit(1);
        }
    }

    // Jeśli po kilku próbach nie udało się osiągnąć dokładności, zakończ program
    if (retries == 0)
    {
        fprintf(stderr, "Nie udało się osiągnąć dokładności po kilku próbach.\n");
        free(resid); free(V); free(workd); free(workl);
        exit(EXIT_FAILURE);
    }

    // Alokacja pamięci dla wyników
    int* select = malloc(ncv * sizeof(int)); // Tablica wyboru wektorów własnych
    memset(select, 0, ncv * sizeof(int));
    check_alloc(select);

    int rvec = 1; // Flaga wskazująca, czy zwracać wektory własne
    float* D = malloc(nev * sizeof(float)); // Wartości własne
    check_alloc(D);
    float sigma = 0.0f; // Przesunięcie spektralne (nieużywane w trybie 1)

    // Wywołanie ARPACK: sseupd_ kończy obliczenia i zwraca wyniki
    sseupd_(&rvec, "A", select, D, V, &n, &sigma, &bmat, &n, which,
            &nev, &tol, resid, &ncv, V, &n, iparam, ipntr,
            workd, workl, &lworkl, &info);

    // Obsługa błędów z sseupd_
    if (info == 1)
    {
        fprintf(stderr, "Błąd ARPACK: Nieosiągnięta dokładność w obliczeniach wartości własnych. Spróbuj zwiększyć ncv lub zmniejszyć tol.\n");
        free(select); free(D); free(V); free(workd); free(workl); free(resid);
        exit(EXIT_FAILURE);
    } else if (info != 0)
    {
        fprintf(stderr, "Błąd ARPACK: %d\n", info);
        free(select); free(D); free(V); free(workd); free(workl); free(resid);
        exit(EXIT_FAILURE);
    }

    // Normalizacja wektorów własnych (każdy wektor jest dzielony przez swoją normę)
    for (int i = 0; i < nev; i++)
    {
        float norm = 0.0f;
        for (int j = 0; j < n; j++)
        {
            float value = V[j * nev + i];
            norm += value * value;
        }
        norm = sqrtf(norm);
        if (norm > 0)
        {
            for (int j = 0; j < n; j++)
            {
                V[j * nev + i] /= norm;
            }
        }
    }

    // Przygotowanie wyników do zwrócenia
    *eigenvectors = malloc(n * nev * sizeof(float)); // Wektory własne
    check_alloc(*eigenvectors);
    *eigenvalues = malloc(nev * sizeof(float)); // Wartości własne
    check_alloc(*eigenvalues);

    // Kopiowanie wyników do tablic wyjściowych
    for (int i = 0; i < nev; i++)
    {
        (*eigenvalues)[i] = D[i];
        for (int j = 0; j < n; j++)
        {
            (*eigenvectors)[i * n + j] = V[j * nev + i];
        }
    }

    // Sortowanie wartości własnych i odpowiadających im wektorów własnych
    for (int i = 0; i < nev - 1; i++)
    {
        for (int j = i + 1; j < nev; j++)
        {
            if ((*eigenvalues)[i] > (*eigenvalues)[j])
            {
                // Zamiana wartości własnych
                float temp_val = (*eigenvalues)[i];
                (*eigenvalues)[i] = (*eigenvalues)[j];
                (*eigenvalues)[j] = temp_val;

                // Zamiana odpowiadających wektorów własnych
                for (int k = 0; k < n; k++)
                {
                    float temp_vec = (*eigenvectors)[i * n + k];
                    (*eigenvectors)[i * n + k] = (*eigenvectors)[j * n + k];
                    (*eigenvectors)[j * n + k] = temp_vec;
                }
            }
        }
    }
    printf("\033[F\033[K");

    // Zwolnienie pamięci
    free(select);
    free(V);
    free(workd);
    free(workl);
    free(D);
    free(resid);
}

// Wypisanie par własnych
void print_eigenpairs(float *eigenvalues, float *eigenvectors, int p, int n)
{
    if(n < max_print_size)
    {
        printf("\n");
        for(int i = 0; i < p; i++)
        {
            printf("\tPara %d: wartość własna = %.6f\n", i, eigenvalues[i]);
            printf("\t  wektor własny[%d] = [", i);
            for(int j = 0; j < n; j++)
            {
                float value = fabs(eigenvectors[i * n + j]) < 1e-6 ? 0.0f : eigenvectors[i * n + j];
                if(j + 1 < n)
                {
                    printf("%.6f, ", value);
                }
                else
                {
                    printf("%.6f", value);
                }
            }
            printf("]\n\n");
        }
    }
    else
    {
        printf("\n\tGraf jest zbyt duży by wypisać pary własne.\n");
    }
}
