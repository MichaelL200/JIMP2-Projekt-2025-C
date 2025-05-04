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
            printf("\tvertex_%d:\n", g_id);
            p_id++;
        }
        else if (current_cluster != clusters[i->vertices_groups[g_id]])
        {
            removed_edges++;
            printf("\t\tUsuwam krawędź: %d -> %d\n", current_vertex, i->vertices_groups[g_id]);

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

    // Usunięcie ostatnich wierzchołków, jeśli są ostatnimi w vertices_ptrs
    while (0 < (int)i->g_count && 0 < (int)i->p_count && i->vertices_ptrs[i->p_count - 1] == (int)i->g_count - 1)
    {
        i->g_count--;
        i->p_count--;
    }

    printf("\n");
    printf("\tUsunięto %d krawędzi\n", removed_edges);

    printf("\n\tGrupy połączeń:\n");
    printv(i->vertices_groups, i->g_count, 10);
    printf("\n\tWskaźniki na pierwsze wierzchołki w grupach:\n");
    printv(i->vertices_ptrs, i->p_count, 10);

    return removed_edges;
}

// Wypisanie rezultatu
void print_result(Result *result)
{
    printf("\n\tRezultat:\n");
    printf("\t\tWynik: %c\n", result->res);
    printf("\t\tPodział: %d\n", result->parts);
    printf("\t\tLiczba usuniętych krawędzi: %d\n", result->cut_count);
    printf("\t\tZachowany margines: %d%%\n", result->margin_kept); // Add percentage symbol
}

// Wypisanie wyników do pliku
void write_output(char *output_file, Input *i, Result *r, int *A, int n, char *format)
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

        // Write max_vertices
        fwrite(&i->max_vertices, sizeof(int), 1, file);

        // Write row_indices
        fwrite(&i->r_count, sizeof(int), 1, file); // Length of row_indices
        fwrite(i->row_indices, sizeof(int), i->r_count, file);

        // Write first_vertices
        fwrite(&i->f_count, sizeof(int), 1, file); // Length of first_vertices
        fwrite(i->first_vertices, sizeof(int), i->f_count, file);

        // Write vertices_groups
        fwrite(&i->g_count, sizeof(size_t), 1, file); // Length of vertices_groups
        fwrite(i->vertices_groups, sizeof(int), i->g_count, file);

        // Write vertices_ptrs
        fwrite(&i->p_count, sizeof(size_t), 1, file); // Length of vertices_ptrs
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