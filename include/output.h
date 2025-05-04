#ifndef OUTPUT_H
#define OUTPUT_H

#include "input.h"

// Struktura do przechowywania wyniku
typedef struct
{
    char res;
    int parts;
    int cut_count;
    int margin_kept;
} Result;

// Wypisanie wyników do pliku
void write_output(char *output_file, Input *i, int *A, int n, char *format);

#endif // OUTPUT_H