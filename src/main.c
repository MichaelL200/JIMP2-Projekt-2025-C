#include <stdio.h>
#include <unistd.h> // dla getopt
#include <stdlib.h> // dla exit

int main(int argc, char *argv[])
{
    int opt;
    int parts = 2;
    int margin = 10;
    char *input_file = NULL;
    char *output_file = "output.txt";
    char *format = "txt";

    // Parsowanie argumentów
    while ((opt = getopt(argc, argv, "p:m:o:f:")) != -1)
    {
        switch (opt)
        {
            // Parts
            case 'p':
                parts = atoi(optarg);
                break;
            // Margin
            case 'm':
                margin = atoi(optarg);
                break;
            // Output
            case 'o':
                output_file = optarg;
                break;
            // Format
            case 'f':
                format = optarg;
                break;
            // Nieznana flaga
            case '?':
                fprintf(stderr, "Nieznana flaga: -%c\n", optopt);
                exit(EXIT_FAILURE);
        }
    }

    // Sprawdzenie, czy plik wejściowy został podany
    if (optind < argc)
    {
        input_file = argv[optind];
    } else
    {
        fprintf(stderr, "Brak pliku wejściowego.\n");
        exit(EXIT_FAILURE);
    }

    // Wyświetlenie wczytanych argumentów
    printf("\tPlik wejściowy: %s\n", input_file);
    printf("\tLiczba części: %d\n", parts);
    printf("\tMargines: %d%%\n", margin);
    printf("\tPlik wyjściowy: %s\n", output_file);
    printf("\tFormat: %s\n", format);

    return 0;
}
