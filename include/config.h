#ifndef CONFIG_H
#define CONFIG_H

// Struktura do przechowywania ustawień
typedef struct
{
    int parts;
    int margin;
    char *input_file;
    char *output_file;
    char *format;
} Config;

// Funkcja sprawdzająca, czy nazwa pliku zawiera tylko dozwolone znaki
int is_valid_filename(char *filename);

// Parsowanie argumentów, zwrócenie skonfigurowanej struktury
Config parse_args(int argc, char **argv);

// Weryfikacja pliku wejściowego – sprawdza czy został podany i czy ma rozszerzenie .csrrg
void validate_input_file(Config *c, int argc, char **argv);

// Wyświetlanie konfiguracji
void print_config(const Config *c);

#endif // CONFIG_H
