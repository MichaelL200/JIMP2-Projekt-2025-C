#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>

#define error(msg) do { \
    fprintf(stderr, "%s\n", msg); \
    exit(EXIT_FAILURE); \
} while (0);

#define check_alloc(ptr) \
    do { \
        if ((ptr) == NULL) { \
            fprintf(stderr, "Błąd: Nie udało się przydzielić pamięci dla zmiennej: %s\n", #ptr); \
            exit(EXIT_FAILURE); \
        } \
    } while (0);

#endif // UTILS_H