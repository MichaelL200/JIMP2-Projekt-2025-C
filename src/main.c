#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"

void printv(int *v, int n, int n_row)
{
    printf("\t[ ");
    for(int i = 0; i < n; i++)
    {
        if( !( i % n_row) )
        {
            printf("\n\t");
        }
        printf("%d ", v[i]);
    }
    printf("]\n");
}

void add_edge(int *v, int a, int b, int n)
{
    // [a, b] = [a * n + b]
    v[a * n + b] = 1;
    v[b * n + a] = 1;
    //symetria względem diagonali
}

int getv(int *v, int a, int b, int n)
{
    return v[a * n + b];
}

int main(int argc, char *argv[])
{
    // Parsowanie argumentów
    Config config = parse_args(argc, argv);

    validate_input_file(&config, argc, argv);

    print_config(&config);

    // Otwarcie pliku wejściowego
    FILE *in = fopen(config.input_file, "r");
    if(in == NULL)
    {
        fprintf(stderr, "\tBłąd otwarcia pliku\n");
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int line_number = 0;
    int max_vertices = 0;
    int *vertices_groups = NULL;
    size_t g_count = 0;
    int *vertices_ptrs = NULL;
    size_t p_count = 0;
    int v_count = 0;

    // Otrzytywanie pliku linia po linii
    while((read = getline(&line, &len, in)) != -1)
    {
        line_number++;

        // Odczytanie maksymalnej liczby wierzchołków w wierszu
        if(line_number == 1)
        {
            max_vertices = atoi(line);
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
                    int *tmp = realloc(vertices_groups, (g_count + 1) * sizeof(int));
                    if(tmp == NULL)
                    {
                        fprintf(stderr, "Błąd alokacji pamięci (vertices_groups)\n");
                        free(vertices_groups);
                        free(vertices_ptrs);
                        free(line);
                        fclose(in);
                        exit(EXIT_FAILURE);
                    }
                    vertices_groups = tmp;
                    vertices_groups[g_count++] = val;
                    token = strtok(NULL, ";");

                    // Przypisz znaleziony większy numer indeksu wierzchołka
                    if(val > v_count)
                    {
                        v_count = val;
                    }
                }
                // Odczytanie wskaźników
                else if(line_number == 5)
                {
                    int *tmp = realloc(vertices_ptrs, (p_count + 1) * sizeof(int));
                    if(tmp == NULL)
                    {
                        fprintf(stderr, "Błąd alokacji pamięci (vertices_groups)\n");
                        free(vertices_groups);
                        free(vertices_ptrs);
                        free(line);
                        fclose(in);
                        exit(EXIT_FAILURE);
                    }
                    vertices_ptrs = tmp;
                    vertices_ptrs[p_count++] = val;
                    token = strtok(NULL, ";");
                }
            }
        }
    }

    // Zwolnienie pamięci
    free(line);

    // Zamknięcie pliku wejściowego
    fclose(in);

    // Wyświetlenie wczytanych danych
    printf("\n\tLimit wierzchołków w wierszu: %d\n", max_vertices);
    printf("\n\tGrupy połączeń:\n");
    printv(vertices_groups, g_count, 10);
    printf("\n\tWskaźniki na pierwsze wierzchołki w grupach:\n");
    printv(vertices_ptrs, p_count, 10);
    printf("\n\tLiczba wierzchołków: %d\n", ++v_count);

    // Wczytanie grafu do macierzy sąsiedztwa A
    printf("\n\tPołączenia dodane do macierzy sąsiedztwa:");
    int *A = calloc(v_count * v_count, sizeof(int));
    int p = 0;
    int v = 0;
    for(int i = 0; i < (int)g_count; i++)
    {
        // Zaktualizuj wierzchołek
        if(p < (int)p_count && i == vertices_ptrs[p])
        {
            v = vertices_groups[i];
            printf("\n\t%d - %d: ", i, v);
            p++;
        }
        else
        {
            printf("\t%d", vertices_groups[i]);
            add_edge(A, v, vertices_groups[i], v_count);
        }
    }
    printf("\n");

    // Test
    //printf("\]n\t%d\n", getv(A, 93, 90, v_count));

    // Wyświetlenie macierzy sąsiedztwa
    printv(A, v_count * v_count, v_count);

    // Zwolnienie pamięci
    free(A);
    free(vertices_ptrs);
    free(vertices_groups);

    return EXIT_SUCCESS;
}
