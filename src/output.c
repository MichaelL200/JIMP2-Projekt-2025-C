#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <string.h>

#include "output.h"
#include "input.h"

// Wypisanie wyników do pliku
void write_output(char *output_file, Input *i, int *A, int n, char *format);