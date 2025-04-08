#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "input.h"
#include "mat_vec.h"

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
        // Pominięcie linii dla interfejsu graficznego
        else if(line_number == 2 || line_number == 3)
        {
            continue;
        }
        // Odczytanie grafu (krawędzi pomiędzy wierzchołkami i ich numerów)
        else if(line_number == 4 || line_number == 5)
        {
            // Usunięcie znaku nowej linii
            if(line[read - 1] == '\n')
            {
                line[read - 1] = '\0';
            }

            // Podzielenie linii na tokeny przy użyciu średnika jako separatora
            char *token = strtok(line,";");
            while(token != NULL)
            {
                int val = atoi(token);

                // Odczytanie grup
                if(line_number == 4)
                {
                    int *tmp = realloc(i->vertices_groups, (i->g_count + 1) * sizeof(int));
                    if(tmp == NULL)
                    {
                        fprintf(stderr, "Błąd alokacji pamięci (vertices_groups)\n");
                        free(line);
                        free(i->vertices_groups);
                        free(i->vertices_ptrs);
                        fclose(i->in);
                        exit(EXIT_FAILURE);
                    }
                    i->vertices_groups = tmp;
                    i->vertices_groups[i->g_count++] = val;
                    token = strtok(NULL, ";");

                    // Przypisz znaleziony większy numer indeksu wierzchołka
                    if(val > i->v_count)
                    {
                        i->v_count = val;
                    }
                }
                // Odczytanie wskaźników
                else if(line_number == 5)
                {
                    int *tmp = realloc(i->vertices_ptrs, (i->p_count + 1) * sizeof(int));
                    if(tmp == NULL)
                    {
                        fprintf(stderr, "Błąd alokacji pamięci (vertices_groups)\n");
                        free(line);
                        free(i->vertices_groups);
                        free(i->vertices_ptrs);
                        fclose(i->in);
                        exit(EXIT_FAILURE);
                    }
                    i->vertices_ptrs = tmp;
                    i->vertices_ptrs[i->p_count++] = val;
                    token = strtok(NULL, ";");
                }
            }
        }
    }

    // Zwolnienie pamięci dla odczytu linii
    free(line);

    // Zamknięcie pliku wejściowego
    fclose(i->in);
}

// Wypisanie informacji wczytanych z pliku wejściowego
void print_input(Input *i)
{
    // Wyświetlenie wczytanych danych
    printf("\n\tLimit wierzchołków w wierszu: %d\n", i->max_vertices);
    printf("\n\tGrupy połączeń:\n");
    printv(i->vertices_groups, i->g_count, 10);
    printf("\n\tWskaźniki na pierwsze wierzchołki w grupach:\n");
    printv(i->vertices_ptrs, i->p_count, 10);
    printf("\n\tLiczba wierzchołków: %d\n", ++(i->v_count));
}