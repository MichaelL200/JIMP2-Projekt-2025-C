#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "input.h"
#include "mat_vec.h"
#include "config.h"
#include "utils.h"

// Otwarcie pliku wejściowego i inicjalizacja
void open_init_input(Input *i, const char *file)
{
    // Inicjalizacja struktury
    memset(i, 0, sizeof(Input));

    // Otwarcie pliku wejściowego
    i->in = fopen(file, "r");
    if(i->in == NULL)
    {
        fprintf(stderr, "\tBłąd otwarcia pliku\n");
        exit(EXIT_FAILURE);
    }
}

// Czytanie pliku linia po linii
void read_input(Input *i)
{
    ssize_t read;
    char *line = NULL;
    int line_number = 0;

    while((read = getline(&line, &i->len, i->in)) != -1)
    {
        line_number++;

        // Odczytanie maksymalnej liczby wierzchołków w wierszu
        if(line_number == 1)
        {
            i->max_vertices = atoi(line);
        }
        // Odczytanie indeksów wierszy
        else if(line_number == 2)
        {
            char *token = strtok(line, ";");
            while(token != NULL)
            {
                int val = atoi(token);
                int *tmp = realloc(i->row_indices, (i->r_count + 1) * sizeof(int));
                if (tmp == NULL)
                check_alloc(tmp);
                i->row_indices = tmp;
                i->row_indices[i->r_count++] = val;
                token = strtok(NULL, ";");
            }
        }
        // Odczytanie pierwszych wierzchołków wierszy
        else if(line_number == 3)
        {
            char *token = strtok(line, ";");
            while(token != NULL)
            {
                int val = atoi(token);
                int *tmp = realloc(i->first_vertices, (i->f_count + 1) * sizeof(int));
                check_alloc(tmp);
                i->first_vertices = tmp;
                i->first_vertices[i->f_count++] = val;
                token = strtok(NULL, ";");
            }
        }
        // Odczytanie grup
        else if(line_number == 4)
        {
            char *token = strtok(line, ";");
            while(token != NULL)
            {
                int val = atoi(token);
                int *tmp = realloc(i->vertices_groups, (i->g_count + 1) * sizeof(int));
                check_alloc(tmp);
                i->vertices_groups = tmp;
                i->vertices_groups[i->g_count++] = val;

                if(val > i->v_count)
                {
                    i->v_count = val;
                }
                token = strtok(NULL, ";");
            }
        }
        // Odczytanie wskaźników
        else if(line_number == 5)
        {
            char *token = strtok(line, ";");
            while(token != NULL)
            {
                int val = atoi(token);
                int *tmp = realloc(i->vertices_ptrs, (i->p_count + 1) * sizeof(int));
                check_alloc(tmp);
                i->vertices_ptrs = tmp;
                i->vertices_ptrs[i->p_count++] = val;
                token = strtok(NULL, ";");
            }
        }
    }

    // Zwolnienie pamięci dla odczytu linii
    free(line);

    // Zamknięcie pliku wejściowego
    fclose(i->in);
}

// Sprawdzenie, czy dane wsadowe są poprawne dla tego grafu
void check_input_data(int parts, int count, int *margin)
{
    // Sprawdzenie, czy liczba podziałów nie jest zbyt duża w porównaniu do liczby wierzchołków
    if (2 * parts > count)
    {
        fprintf(stderr, "Błąd: Podwojona liczba podziałów (2 * %d > %d) jest większa niż liczba wierzchołków (%d). Metoda spektralna nieoptymalna.\n", parts, 2 * parts, count);
        exit(EXIT_FAILURE);
    }

    // Obliczenie minimalnego marginesu
    int min_margin = count / parts - (count + parts - 1) / parts;

    // Zwiększenie minimalnego marginesu o 10% lub o stałą wartość
    int adjusted_margin = min_margin + (min_margin / 10); // Dodaj 10% minimalnego marginesu
    if (adjusted_margin == min_margin)
    {
        adjusted_margin += 10; // Jeśli 10% to 0, zwiększ margines o 10
    }

    // Sprawdzenie, czy margines jest wystarczająco duży
    if (*margin < adjusted_margin)
    {
        fprintf(stderr, "Błąd: Margines (%d) jest zbyt mały. Minimalny margines dla %d części i %d wierzchołków wynosi %d.\n", *margin, parts, count, adjusted_margin);
        *margin = adjusted_margin; // Ustaw margines na zwiększony minimalny
        fprintf(stderr, "Ustawiono margines na minimalny powiększony: %d.\n", *margin);
    }
}

// Wypisanie informacji wczytanych z pliku wejściowego
void print_input(Input *i)
{
    // Wyświetlenie wczytanych danych
    printf("\n\tLimit wierzchołków w wierszu: %d\n", i->max_vertices);

    if (i->v_count >  max_print_size) {
        printf("\n\tGraf jest zbyt duży, by wyświetlić wczytane dane grafu.\n");
    } else {
        printf("\n\tIndeksy wierszy:\n");
        printv(i->row_indices, i->r_count, 10);
        printf("\n\tPierwsze wierzchołki w wierszach:\n");
        printv(i->first_vertices, i->f_count, 10);
        printf("\n\tGrupy połączeń:\n");
        printv(i->vertices_groups, i->g_count, 10);
        printf("\n\tWskaźniki na pierwsze wierzchołki w grupach:\n");
        printv(i->vertices_ptrs, i->p_count, 10);
    }

    printf("\n\tLiczba wierzchołków: %d\n", ++(i->v_count));
}

// Funkcja zwalniająca pamięć dla struktury Input
void free_input(Input *i)
{
    if(i == NULL)
    {
        error("Zmienna Input nie może być NULL");
    }

    if (i->row_indices) free(i->row_indices);
    if (i->first_vertices) free(i->first_vertices);
    if (i->vertices_groups) free(i->vertices_groups);
    if (i->vertices_ptrs) free(i->vertices_ptrs);
}
