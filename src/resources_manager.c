#include <stdlib.h>
#include <stdio.h>

#include "resource_manager.h"

// Deklaracja zewnętrznego wskaźnika do kontekstu zasobów
ResourceContext *res_ctx = NULL;

// Inicjalizacja kontekstu zasobów
void init_resource_context()
{
    res_ctx = (ResourceContext *)malloc(sizeof(ResourceContext));
    if (res_ctx == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla kontekstu zasobów!\n");
        exit(EXIT_FAILURE);
    }
    res_ctx->capacity = 5;
    res_ctx->count = 0;
    res_ctx->resources = (void **)malloc(res_ctx->capacity * sizeof(void *));
    if (res_ctx->resources == NULL)
    {
        fprintf(stderr, "Błąd alokacji pamięci dla kontekstu zasobów!\n");
        exit(EXIT_FAILURE);
    }
}

// Zwalnianie wszystkich zarejestrowanych zasobów
void cleanup_resources()
{
    for (size_t i = 0; i < res_ctx->count; i++)
    {
        if(res_ctx->resources[i] != NULL)
        {
            free(res_ctx->resources[i]);
        }
    }
    free(res_ctx->resources);
    res_ctx->resources = NULL;
    res_ctx->count = 0;
    res_ctx->capacity = 0;
    free(res_ctx);
}

// Rejestrowanie zasobów (wskaźników)
void register_resource(void *ptr)
{
    if (res_ctx == NULL)
    {
        fprintf(stderr, "Błąd: nieprawidłowy wskaźnik kontekstu zasobów!\n");
        cleanup_resources();
        exit(EXIT_FAILURE);
    }
    if(ptr == NULL)
    {
        fprintf(stderr, "Błąd: nieprawidłowy wskaźnik zasobu!\n");
        cleanup_resources();
        exit(EXIT_FAILURE);
    }

    // Jeśli brak miejsca w tablicy resources
    if (res_ctx->count == res_ctx->capacity)
    {
        // Zwiększenie pojemności tablicy resources
        size_t new_capacity = res_ctx->capacity * 2;
        void **new_resources = (void **)realloc(res_ctx->resources, new_capacity * sizeof(void *));
        if (new_resources == NULL)
        {
            fprintf(stderr, "Błąd alokacji pamięci dla zasobów!\n(zwiększenia pojemności tablicy\n");
            cleanup_resources();
            exit(EXIT_FAILURE);
        }
        res_ctx->resources = new_resources;
        res_ctx->capacity = new_capacity;
    }

    // Dodanie wskaźnika do tablicy resources
    res_ctx->resources[res_ctx->count++] = ptr;
}

// Usuwanie zasobu z kontekstu
void deregister_resource(void *ptr)
{
    if (res_ctx == NULL)
    {
        fprintf(stderr, "Błąd: nieprawidłowy wskaźnik kontekstu zasobów!\n");
        cleanup_resources();
        exit(EXIT_FAILURE);
    }
    if(ptr == NULL)
    {
        fprintf(stderr, "Błąd: nieprawidłowy wskaźnik zasobu!\n");
        cleanup_resources();
        exit(EXIT_FAILURE);
    }

    // Szukanie wskaźnika w tablicy resources
    size_t i;
    for (i = 0; i < res_ctx->count; i++)
    {
        if (res_ctx->resources[i] == ptr)
        {
            // Zwalnianie pamięci zasobu
            free(res_ctx->resources[i]);

            // Przesuwanie pozostałych wskaźników w tablicy
            for (size_t j = i; j < res_ctx->count - 1; j++)
            {
                res_ctx->resources[j] = res_ctx->resources[j + 1];
            }

            // Zmniejszenie licznika zasobów
            res_ctx->count--;
            return;
        }
    }

    // Jeśli zasób nie został znaleziony
    fprintf(stderr, "Błąd: zasób nie znaleziony w kontekście!\n");
    cleanup_resources();
    exit(EXIT_FAILURE);
}
