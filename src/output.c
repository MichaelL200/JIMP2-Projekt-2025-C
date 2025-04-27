#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <string.h>

#include "output.h"
#include "input.h"

int *resize_array(int *array, int current_capacity, int new_capacity)
{
    int *new_array = realloc(array, new_capacity * sizeof(int));
    if (new_array == NULL)
    {
        fprintf(stderr, "Błąd: Nie można rozszerzyć pamięci dla tablicy.\n");
        free(array);
        exit(EXIT_FAILURE);
    }

    // Inicjalizacja nowo przydzielonej pamięci
    memset(new_array + current_capacity, 0, (new_capacity - current_capacity) * sizeof(int));

    return new_array;
}

// Wypisywanie grafu do pliku
void write_output(char *output_file, Result *r, Input *i, int *A, int n, char *format)
{
    // Sprawdzenie, czy format jest poprawny
    if (strcmp(format, "txt") != 0 && strcmp(format, "bin") != 0)
    {
        fprintf(stderr, "Niepoprawny format pliku wyjściowego. Użyj 'txt' lub 'bin'.\n");
        exit(EXIT_FAILURE);
    }

    // Tworzenie folderu output, jeśli nie istnieje
    if (mkdir("output", 0777) == -1 && errno != EEXIST)
    {
        fprintf(stderr, "Nie można utworzyć folderu output\n");
        exit(EXIT_FAILURE);
    }

    // Tworzenie ścieżki do pliku wyjściowego w folderze output/
    char output_path[256];
    snprintf(output_path, sizeof(output_path), "output/%s", output_file);

    // Otworzenie pliku wyjściowego
    FILE *out = fopen(output_path, "w");
    if (out == NULL)
    {
        fprintf(stderr, "Nie można otworzyć pliku wyjściowego: %s\n", output_path);
        exit(EXIT_FAILURE);
    }

    if(strcmp(format, "txt") == 0)
    {
        // Wypisanie wyniku działania programu do pliku (0)
        fprintf(out, "%c ", r->res);
        fprintf(out, "%d ", r->parts);
        fprintf(out, "%d ", r->cut_count);
        fprintf(out, "%d\n", r->margin_kept);

        // Wypisanie maksymalnej liczby wierzchołków w wierszu (1)
        fprintf(out, "%d\n", i->max_vertices);

        // Wypisanie indeksów wierszy (2)
        for (int j = 0; j < i->r_count; j++)
        {
            fprintf(out, "%d", i->row_indices[j]);
            if (j < i->r_count - 1)
            {
                fprintf(out, ";");
            }
        }
        fprintf(out, "\n");

        // Wypisanie pierwszych wierzchołków w wierszach (3)
        for (int j = 0; j < i->f_count; j++)
        {
            fprintf(out, "%d", i->first_vertices[j]);
            if (j < i->f_count - 1)
            {
                fprintf(out, ";");
            }
        }
        fprintf(out, "\n");

        // Wypisanie grup wierzchołków z A (4) i ich przypisania na podstawie wskaźników (5)
        int index = 0;
        int count = 0;
        int *indices = (int *)malloc(n * sizeof(int));
        if (indices == NULL)
        {
            fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla wskaźników.\n");
            fclose(out);
            exit(EXIT_FAILURE);
        }
        // Wypisanie grup wierzchołków z A (4)
        for (int r = 0; r < n; r++)
        {
            indices[count++] = index++;
            fprintf(out, "%d;", r); // Wypisanie wierzchołka
            for (int c = r + 1; c < n; c++) // Tylko elementy po prawej stronie diagonali
            {
                if (A[r * n + c] == 1)
                {
                    fprintf(out, "%d;", c); // Wypisanie wierzchołka
                    index++;
                }
            }
        }
        fprintf(out, "\n");
        // Wypisanie wskaźników (5)
        for (int j = 0; j < n; j++)
        {
            fprintf(out, "%d", indices[j]);
            if (j < n - 1)
            {
                fprintf(out, ";");
            }
        }
    }
    else if (strcmp(format, "bin") == 0)
    {
        // Wypisanie wyniku działania programu do pliku (0)
        fwrite(&r->res, sizeof(char), 1, out);
        fwrite(&r->parts, sizeof(int), 1, out);
        fwrite(&r->cut_count, sizeof(int), 1, out);
        fwrite(&r->margin_kept, sizeof(int), 1, out);

        // Wypisanie maksymalnej liczby wierzchołków w wierszu (1)
        fwrite(&i->max_vertices, sizeof(int), 1, out);

        // Wypisanie indeksów wierszy (2)
        fwrite(i->row_indices, sizeof(int), i->r_count, out);

        // Wypisanie pierwszych wierzchołków w wierszach (3)
        fwrite(i->first_vertices, sizeof(int), i->f_count, out);

        // Zapisanie grup wierzchołków z A (4) i ich przypisania na podstawie wskaźników (5)
        int index = 0;
        int count = 0;
        int capacity = n;
        int *groups = (int *)malloc(capacity * sizeof(int));
        if (groups == NULL)
        {
            fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla grup wierzchołków.\n");
            free(groups);
            fclose(out);
            exit(EXIT_FAILURE);
        }
        int *indices = (int *)malloc(n * sizeof(int));
        if (indices == NULL)
        {
            fprintf(stderr, "Błąd: Nie można przydzielić pamięci dla wskaźników.\n");
            free(groups);
            free(indices);
            fclose(out);
            exit(EXIT_FAILURE);
        }
        // Zapisanie grup wierzchołków z A (4)
        for (int r = 0; r < n; r++)
        {
            if (index >= capacity)
            {
                capacity *= 2;
                groups = resize_array(groups, capacity / 2, capacity);
            }
            indices[count++] = index;
            groups[index++] = r;
            for (int c = r + 1; c < n; c++)
            {
                if (A[r * n + c] == 1)
                {
                    if (index >= capacity)
                    {
                        capacity *= 2;
                        groups = resize_array(groups, capacity / 2, capacity);
                    }
                    groups[index++] = c;
                }
            }
        }

        // Wypisanie grup wierzchołków z A (4)
        fwrite(&index, sizeof(int), 1, out); // Zapisz liczbę elementów w groups
        fwrite(groups, sizeof(int), index, out);

        // Wypisanie wskaźników (5)
        fwrite(&n, sizeof(int), 1, out); // Zapisz liczbę elementów w indices
        fwrite(indices, sizeof(int), n, out);

        free(groups);
        free(indices);
    }

    fclose(out);
}