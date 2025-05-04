#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <string.h>

#include "output.h"
#include "input.h"
#include "mat_vec.h"

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

            g_id--; // Cofanie indeksu, by nie pominąć następnego elementu
        }
        //printf("\ti->vertices_groups[%d] = %d\n", g_id, i->vertices_groups[g_id]);
    }

    // Usunięcie ostatnich wierzchołków, jeśli są ostatnimi w vertices_ptrs
    while (i->g_count > 0 && i->p_count > 0 && i->vertices_ptrs[i->p_count - 1] == i->g_count - 1)
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

// Wypisanie wyników do pliku
void write_output(char *output_file, Input *i, int *A, int n, char *format);
