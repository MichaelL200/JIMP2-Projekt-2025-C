#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>

// Struktura do przechowywania danych z pliku wejściowego
typedef struct
{
    FILE *in;
    int max_vertices;
    int *vertices_groups;
    size_t g_count;
    int *vertices_ptrs;
    size_t p_count;
    int v_count;
    size_t len;
} Input;

// Otwarcie pliku wejściowego i inicjalizacja
void open_init_input(Input *i, const char *file);

// Czytanie pliku linia po linii
void read_input(Input *i);

// Wypisanie informacji wczytanych z pliku wejściowego
void print_input(Input *i);

#endif // INPUT_H
