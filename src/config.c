#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "config.h"

// Funkcja sprawdzająca, czy nazwa pliku zawiera tylko dozwolone znaki
int is_valid_filename(char *filename)
{
    // Sprawdzenie, czy nazwa pliku zawiera tylko dozwolone znaki
    for (int i = 0; filename[i] != '\0'; i++)
    {
        if (filename[i] == '/' || filename[i] == '\\' || filename[i] == ':' || filename[i] == '*' || filename[i] == '?' || filename[i] == '"' || filename[i] == '<' || filename[i] == '>' || filename[i] == '|')
        {
            return 0; // Niedozwolony znak
        }
    }
    return 1; // Dozwolone znaki
}

// Parsowanie argumentów, zwrócenie skonfigurowanej struktury
Config parse_args(int argc, char **argv)
{
    Config c;

    // Ustawienia domyślne
    c.parts = 2;
    c.margin = 10;
    c.input_file = NULL;
    c.output_file = "output.txt";
    c.format = "txt";

    // Parsowanie argumentów
    int opt;
    while ((opt = getopt(argc, argv, "p:m:o:f:")) != -1)
    {
        switch (opt)
        {
            // Liczba części
            case 'p':
                c.parts = atoi(optarg);
                if (c.parts < 2)
                {
                    fprintf(stderr, "Liczba części musi być większa od 2\n");
                    c.parts = 2;
                }
                break;
            // Margines
            case 'm':
                c.margin = atoi(optarg);
                if (c.margin < 0)
                {
                    fprintf(stderr, "Margines musi być większy od 0.\n");
                    c.margin = 10;
                }
                break;
            // Plik wyjściowy
            case 'o':
                if (is_valid_filename(optarg))
                {
                    c.output_file = optarg;
                }
                else
                {
                    fprintf(stderr, "Nazwa pliku wyjściowego zawiera niedozwolone znaki. Ustawiono domyślną nazwę: output.txt\n");
                    c.output_file = "output.txt";
                }
                break;
            // Format
            case 'f':
                c.format = optarg;
                // Sprawdzenie poprawności formatu
                if (strcmp(c.format, "txt") != 0 && strcmp(c.format, "bin") != 0)
                {
                    fprintf(stderr, "Nieznany format: %s. Ustawiono domyślny format: txt\n", c.format);
                    c.format = "txt";
                }
                break;
            // Nieznana flaga
            case '?':
                fprintf(stderr, "Nieznana flaga: -%c\n", optopt);
        }
    }

    // Sprawdzenie rozszerzenia pliku wyjściowego i poprawienie go, jeśli jest nieprawidłowe
    if (c.output_file != NULL)
    {
        const char *dot = strrchr(c.output_file, '.'); // Znajdź ostatnią kropkę w nazwie pliku
        if (dot == NULL || strcmp(dot + 1, c.format) != 0)
        {
            // Nowe miejsce na nazwę z poprawnym rozszerzeniem
            size_t base_len = dot ? (size_t)(dot - c.output_file) : strlen(c.output_file);
            size_t ext_len = strlen(c.format);
            char *new_output = malloc(base_len + ext_len + 2); // +1 na '.' +1 na '\0'
            if (!new_output)
            {
                perror("malloc");
                exit(EXIT_FAILURE);
            }
            strncpy(new_output, c.output_file, base_len); // Skopiuj nazwę pliku bez rozszerzenia
            new_output[base_len] = '.';                  // Dodaj kropkę
            strcpy(new_output + base_len + 1, c.format); // Dodaj nowe rozszerzenie
            c.output_file = new_output;                 // Zaktualizuj nazwę pliku wyjściowego
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
        if(len_f < len_e || strcmp(c->input_file + len_f - len_e, ext) != 0)
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
