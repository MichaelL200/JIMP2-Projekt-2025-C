#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>

// Definicja maksymalnego rozmiaru grafu do wypisywania danych
extern int max_print_size;

// Struktura do przechowywania danych z pliku wejściowego
typedef struct
{
    FILE *in;
    int max_vertices;
    int *row_indices;
    int r_count;
    int *first_vertices;
    int f_count;
    int *vertices_groups;
    size_t g_count;
    int *vertices_ptrs;
    size_t p_count;
    int v_count;
    size_t len;
} Input;

// Funkcja sprawdzająca, czy nazwa pliku zawiera tylko dozwolone znaki
int is_valid_filename(char *filename);

// Otwarcie pliku wejściowego i inicjalizacja
void open_init_input(Input *i, const char *file);

// Czytanie pliku linia po linii
void read_input(Input *i);

// Sprawdzenie, czy dane wsadowe są poprawne dla tego grafu
void check_input_data(int parts, int count, int *margin);

// Wypisanie informacji wczytanych z pliku wejściowego
void print_input(Input *i);

// Funkcja zwalniająca pamięć dla struktury Input
void free_input(Input *i);

#endif // INPUT_H
