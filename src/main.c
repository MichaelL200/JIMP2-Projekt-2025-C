#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

void printv(int *v, int n)
{
    printf("\n\t[ ");
    for(int i = 0; i < n; i++)
    {
        if( !( (i + 1) % 10) )
        {
            printf("\n\t");
        }
        printf("%d ", v[i]);
    }
    printf("]\n");
}

int main(int argc, char *argv[])
{
    int opt;
    int parts = 2;
    int margin = 10;
    char *input_file = NULL;
    char *output_file = "output.txt";
    char *format = "txt";

    // Parsowanie argumentów
    while ((opt = getopt(argc, argv, "p:m:o:f:")) != -1)
    {
        switch (opt)
        {
            // Parts
            case 'p':
                parts = atoi(optarg);
                break;
            // Margin
            case 'm':
                margin = atoi(optarg);
                break;
            // Output
            case 'o':
                output_file = optarg;
                break;
            // Format
            case 'f':
                format = optarg;
                break;
            // Nieznana flaga
            case '?':
                fprintf(stderr, "\tNieznana flaga: -%c\n", optopt);
                exit(EXIT_FAILURE);
        }
    }

    // Sprawdzenie, czy plik wejściowy został podany
    if (optind < argc)
    {
        input_file = argv[optind];

        // Sprawdzenie rozszerzenia pliku wejściowego
        const char* ext = ".csrrg";
        size_t len_f = strlen(input_file);
        size_t len_e = strlen(ext);
        if(len_f < len_e ||
         strcmp(input_file + len_f - len_e, ext) != 0)
        {
            fprintf(stderr, "\tZłe rozszerzenie pliku wejściowego. Popawne rozszerzenie to .csrrg\n");
            exit(EXIT_FAILURE);
        }
        
    } else
    {
        fprintf(stderr, "\tBrak pliku wejściowego\n");
        exit(EXIT_FAILURE);
    }

    // Wyświetlenie wczytanych argumentów
    printf("\tPlik wejściowy: %s\n", input_file);
    printf("\tLiczba części: %d\n", parts);
    printf("\tMargines: %d%%\n", margin);
    printf("\tPlik wyjściowy: %s\n", output_file);
    printf("\tFormat: %s\n", format);

    // Otwarcie pliku wejściowego
    FILE *in = fopen(input_file, "r");
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

    // Zamknięcie pliku wejściowego
    fclose(in);

    // Wyświetlenie wczytanych danych
    printf("\tLimit wierzchołków w wierszu: %d\n", max_vertices);
    printv(vertices_groups, g_count);
    printv(vertices_ptrs, p_count);

    return EXIT_SUCCESS;
}
