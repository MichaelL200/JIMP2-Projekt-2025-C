#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "output.h"
#include "input.h"
#include "mat_vec.h"
#include "utils.h"

// Modyfikacja grafu w oparciu o podział
int modify_graph(Input *i, int* clusters)
{
    int removed_edges = 0;

    printf("\n");
    int p_id = 0;
    int current_vertex;
    int current_cluster;
    for (int g_id = 0; g_id < (int)i->g_count; g_id++)
    {
        if (p_id < (int)i->p_count && i->vertices_ptrs[p_id] == g_id)
        {
            current_vertex = i->vertices_groups[g_id];
            current_cluster = clusters[current_vertex];
            if(i->v_count < max_print_size)
            {
                printf("\tvertex_%d:\n", g_id);
            }
            p_id++;
        }
        else if (current_cluster != clusters[i->vertices_groups[g_id]])
        {
            removed_edges++;
            if(i->v_count < max_print_size)
            {
                printf("\t\tUsuwam krawędź: %d -> %d\n", current_vertex, i->vertices_groups[g_id]);
            }
            // Usunięcie krawędzi z vertices_groups
            for (int j = g_id; j < (int)i->g_count - 1; j++)
            {
                i->vertices_groups[j] = i->vertices_groups[j + 1];
            }
            i->g_count--;

            // Dopasowanie wskaźników vertices_ptrs
            for (int j = p_id; j < (int)i->p_count; j++)
            {
                i->vertices_ptrs[j]--;
            }

            g_id--; // Cofnięcie indeksu, aby nie pominąć następnego elementu
        }
    }

    if(i->v_count > max_print_size)
    {
        printf("\tGraf jest zbyt duży by wypisać wierzchołki i usunięte krawędzie.\n");
    }


    // Usunięcie ostatnich wierzchołków, jeśli są ostatnimi w vertices_ptrs
    while (0 < (int)i->g_count && 0 < (int)i->p_count && i->vertices_ptrs[i->p_count - 1] == (int)i->g_count - 1)
    {
        i->g_count--;
        i->p_count--;
    }

    printf("\n");
    printf("\tUsunięto %d krawędzi\n", removed_edges);

    if (i->v_count < max_print_size)
    {
        printf("\n\tGrupy połączeń:\n");
        printv(i->vertices_groups, i->g_count, 10);
        printf("\n\tWskaźniki na pierwsze wierzchołki w grupach:\n");
        printv(i->vertices_ptrs, i->p_count, 10);
    }
    else
    {
        printf("\n\tGraf jest zbyt duży by wypisać grupy połączeń i wskaźniki po usunięciu krawędzi.\n");
    }

    return removed_edges;
}

// Wypisanie rezultatu
void print_result(Result *result)
{
    printf("\n\tRezultat:\n");
    printf("\t\tWynik: %c\n", result->res);
    printf("\t\tPodział: %d\n", result->parts);
    printf("\t\tLiczba usuniętych krawędzi: %d\n", result->cut_count);
    printf("\t\tZachowany margines: %d%%\n", result->margin_kept); // Zawsze wyświetl margines
    if (result->res == 'F') {
        printf("\t\t(Uwaga: Margines został przekroczony)\n"); // Dodaj ostrzeżenie, jeśli margines przekroczony
    }
}

// Wypisanie wyników do pliku
void write_output(char *output_file, Input *i, Result *r, char *format)
{
    // Upewnij się, że katalog output/ istnieje
    struct stat st = {0};
    if (stat("output", &st) == -1)
    {
        if (mkdir("output", 0700) != 0)
        {
            error("Nie można utworzyć katalogu output/");
        }
    }

    // Dodanie ścieżki katalogu output/ do nazwy pliku
    char full_path[256];
    snprintf(full_path, sizeof(full_path), "output/%s", output_file);

    FILE *file;
    if (strcmp(format, "txt") == 0)
    {
        file = fopen(full_path, "w");
        if (!file)
        {
            error("Nie można otworzyć pliku do zapisu");
        }
        fprintf(file, "%c %d %d %d\n", r->res, r->parts, r->cut_count, r->margin_kept);

        // Wypisanie maksymalnej liczby wierzchołków w wierszu
        fprintf(file, "%d\n", i->max_vertices);

        // Wypisanie indeksów wierszy
        for (int j = 0; j < i->r_count; j++)
        {
            fprintf(file, "%d", i->row_indices[j]);
            if (j < i->r_count - 1)
            {
                fprintf(file, ";");
            }
        }
        fprintf(file, "\n");

        // Wypisanie pierwszych wierzchołków wierszy
        for (int j = 0; j < i->f_count; j++)
        {
            fprintf(file, "%d", i->first_vertices[j]);
            if (j < i->f_count - 1)
            {
                fprintf(file, ";");
            }
        }
        fprintf(file, "\n");

        // Wypisanie grup wierzchołków
        for (size_t j = 0; j < i->g_count; j++)
        {
            fprintf(file, "%d", i->vertices_groups[j]);
            if (j < i->g_count - 1)
            {
                fprintf(file, ";");
            }
        }
        fprintf(file, "\n");

        // Wypisanie wskaźników na pierwsze wierzchołki w grupach
        for (size_t j = 0; j < i->p_count; j++)
        {
            fprintf(file, "%d", i->vertices_ptrs[j]);
            if (j < i->p_count - 1)
            {
                fprintf(file, ";");
            }
        }
        fprintf(file, "\n");
    }
    else if (strcmp(format, "bin") == 0)
    {
        file = fopen(full_path, "wb");
        if (!file)
        {
            error("Nie można otworzyć pliku do zapisu");
        }
        fwrite(&r->res, sizeof(char), 1, file);
        fwrite(&r->parts, sizeof(int), 1, file);
        fwrite(&r->cut_count, sizeof(int), 1, file);
        fwrite(&r->margin_kept, sizeof(int), 1, file);

        // Zapisz max_vertices
        fwrite(&i->max_vertices, sizeof(int), 1, file);

        // Zapisz row_indices
        fwrite(&i->r_count, sizeof(int), 1, file); // Długość row_indices
        fwrite(i->row_indices, sizeof(int), i->r_count, file);

        // Zapisz first_vertices
        fwrite(&i->f_count, sizeof(int), 1, file); // Długość first_vertices
        fwrite(i->first_vertices, sizeof(int), i->f_count, file);

        // Zapisz vertices_groups
        fwrite(&i->g_count, sizeof(size_t), 1, file); // Długość vertices_groups
        fwrite(i->vertices_groups, sizeof(int), i->g_count, file);

        // Zapisz vertices_ptrs
        fwrite(&i->p_count, sizeof(size_t), 1, file); // Długość vertices_ptrs
        fwrite(i->vertices_ptrs, sizeof(int), i->p_count, file);
    }
    else
    {
        fprintf(stderr, "Nieobsługiwany format pliku: %s\n", format);
        exit(EXIT_FAILURE);
    }

    printf("\n\tWyniki zapisano do pliku: %s\n", full_path);

    fclose(file);
}

// Wczytanie wyników z pliku wyjściowego (test)
void bin_read(const char *filename, int n)
{
    // Utwórz lokalne zmienne
    Result result;
    Input input = {NULL, 0, NULL, 0, NULL, 0, NULL, 0, NULL, 0, 0, 0};

    const char *ext = strrchr(filename, '.');
    if (!ext || strcmp(ext, ".bin") != 0)
    {
        return;
    }

    FILE *file = fopen(filename, "rb");
    if (!file)
    {
        fprintf(stderr, "Nie można otworzyć pliku do odczytu: %s\n", filename);
        perror("Szczegóły błędu");
        return;
    }

    // Odczytaj dane z pliku
    fread(&result.res, sizeof(char), 1, file);
    fread(&result.parts, sizeof(int), 1, file);
    fread(&result.cut_count, sizeof(int), 1, file);
    fread(&result.margin_kept, sizeof(int), 1, file);

    fread(&input.max_vertices, sizeof(int), 1, file);

    fread(&input.r_count, sizeof(int), 1, file);
    input.row_indices = malloc(input.r_count * sizeof(int));
    check_alloc(input.row_indices);
    fread(input.row_indices, sizeof(int), input.r_count, file);

    fread(&input.f_count, sizeof(int), 1, file);
    input.first_vertices = malloc(input.f_count * sizeof(int));
    check_alloc(input.first_vertices);
    fread(input.first_vertices, sizeof(int), input.f_count, file);

    fread(&input.g_count, sizeof(size_t), 1, file);
    input.vertices_groups = malloc(input.g_count * sizeof(int));
    check_alloc(input.vertices_groups);
    fread(input.vertices_groups, sizeof(int), input.g_count, file);

    fread(&input.p_count, sizeof(size_t), 1, file);
    input.vertices_ptrs = malloc(input.p_count * sizeof(int));
    check_alloc(input.vertices_ptrs);
    fread(input.vertices_ptrs, sizeof(int), input.p_count, file);

    fclose(file);

    // Wyświetl dane wejściowe
    printf("\n\tDane zapisane w formacie binarnym:\n");
    printf("\t\tWynik: %c\n", result.res);
    printf("\t\tPodział: %d\n", result.parts);
    printf("\t\tLiczba usuniętych krawędzi: %d\n", result.cut_count);
    printf("\t\tZachowany margines: %d%%\n", result.margin_kept);

    printf("\tLimit wierzchołków w wierszu: %d\n", input.max_vertices);
    if(n < max_print_size)
    {
        printf("\tIndeksy wierszy:\n");
        printv(input.row_indices, input.r_count, 10);
        printf("\tPierwsze wierzchołki w wierszach:\n");
        printv(input.first_vertices, input.f_count, 10);
        printf("\tGrupy połączeń:\n");
        printv(input.vertices_groups, input.g_count, 10);
        printf("\tWskaźniki na pierwsze wierzchołki w grupach:\n");
        printv(input.vertices_ptrs, input.p_count, 10);
    }
    else
    {
        printf("\tGraf jest zbyt duży by wypisać tablice z pliku binarnego.\n");
    }

    // Zwolnij zaalokowaną pamięć
    free(input.row_indices);
    free(input.first_vertices);
    free(input.vertices_groups);
    free(input.vertices_ptrs);
}