#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

#include "output.h"
#include "input.h"

// Wypisywanie grafu do pliku
void write_output(char *output_file, Input *i, int *A, int n)
{
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
        perror("Error opening output file");
        exit(EXIT_FAILURE);
    }

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