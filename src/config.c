#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "config.h"

// Parsowanie argumentów, zwrócenie skonfigurowanej struktury
Config parse_args(int argc, char **argv) {
    Config config;

    // Ustawienia domyślne
    config.parts = 2;
    config.margin = 10;
    config.input_file = NULL;
    config.output_file = "output.txt";
    config.format = "txt";

    // Parsowanie argumentów
    int opt;
    while ((opt = getopt(argc, argv, "p:m:o:f:")) != -1) {
        switch (opt) {
            // Parts
            case 'p':
                config.parts = atoi(optarg);
                break;
            // Margin
            case 'm':
                config.margin = atoi(optarg);
                break;
            // Output
            case 'o':
                config.output_file = optarg;
                break;
            // Format
            case 'f':
                config.format = optarg;
                break;
            // Nieznana flaga
            case '?':
                fprintf(stderr, "Nieznana flaga: -%c\n", optopt);
                exit(EXIT_FAILURE);
        }
    }
    return config;
}

// Weryfikacja pliku wejściowego – sprawdza czy został podany i czy ma rozszerzenie .csrrg
void validate_input_file(Config *config, int argc, char **argv)
{
    // Sprawdzenie, czy plik wejściowy został podany
    if (optind < argc)
    {
        config->input_file = argv[optind];

        // Sprawdzenie rozszerzenia pliku wejściowego
        const char* ext = ".csrrg";
        size_t len_f = strlen(config->input_file);
        size_t len_e = strlen(ext);
        if(len_f < len_e ||
         strcmp(config->input_file + len_f - len_e, ext) != 0)
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
void print_config(const Config *config)
{
    // Wyświetlenie wczytanych argumentów
    printf("\tPlik wejściowy: %s\n", config->input_file);
    printf("\tLiczba części: %d\n", config->parts);
    printf("\tMargines: %d%%\n", config->margin);
    printf("\tPlik wyjściowy: %s\n", config->output_file);
    printf("\tFormat: %s\n", config->format);
}
