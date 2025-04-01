#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "config.h"

// Parsowanie argumentów, zwrócenie skonfigurowanej struktury
Config parse_args(int argc, char **argv) {
    Config c;

    // Ustawienia domyślne
    c.parts = 2;
    c.margin = 10;
    c.input_file = NULL;
    c.output_file = "output.txt";
    c.format = "txt";

    // Parsowanie argumentów
    int opt;
    while ((opt = getopt(argc, argv, "p:m:o:f:")) != -1) {
        switch (opt) {
            // Parts
            case 'p':
                c.parts = atoi(optarg);
                break;
            // Margin
            case 'm':
                c.margin = atoi(optarg);
                break;
            // Output
            case 'o':
                c.output_file = optarg;
                break;
            // Format
            case 'f':
                c.format = optarg;
                break;
            // Nieznana flaga
            case '?':
                fprintf(stderr, "Nieznana flaga: -%c\n", optopt);
                exit(EXIT_FAILURE);
        }
    }
    return c;
}

// Weryfikacja pliku wejściowego – sprawdza czy został podany i czy ma rozszerzenie .csrrg
void validate_input_file(Config *c, int argc, char **argv)
{
    // Sprawdzenie, czy plik wejściowy został podany
    if (optind < argc)
    {
        c->input_file = argv[optind];

        // Sprawdzenie rozszerzenia pliku wejściowego
        const char* ext = ".csrrg";
        size_t len_f = strlen(c->input_file);
        size_t len_e = strlen(ext);
        if(len_f < len_e ||
         strcmp(c->input_file + len_f - len_e, ext) != 0)
        {
            fprintf(stderr, "\tZłe rozszerzenie pliku wejściowego. Popawne rozszerzenie to .csrrg\n");
            exit(EXIT_FAILURE);
        }
        
    }
    else
    {
        fprintf(stderr, "\tBrak pliku wejściowego\n");
        exit(EXIT_FAILURE);
    }
}

// Wyświetlanie konfiguracji
void print_config(const Config *c)
{
    // Wyświetlenie wczytanych argumentów
    printf("\tPlik wejściowy: %s\n", c->input_file);
    printf("\tLiczba części: %d\n", c->parts);
    printf("\tMargines: %d%%\n", c->margin);
    printf("\tPlik wyjściowy: %s\n", c->output_file);
    printf("\tFormat: %s\n", c->format);
}
