#ifndef RESOURCE_MANAGER_H
#define RESOURCE_MANAGER_H

#include <stdlib.h>

// Struktura do przechowywania kontekstu zasobów
typedef struct
{
    void **resources;
    size_t count;
    size_t capacity; 
} ResourceContext;

// Deklaracja zewnętrznego wskaźnika do kontekstu zasobów
extern ResourceContext *res_ctx;

// Inicjalizacja kontekstu zasobów
void init_resource_context();

// Zwalnianie wszystkich zarejestrowanych zasobów
void cleanup_resources();

// Rejestrowanie zasobów (wskaźników)
void register_resource(void *ptr);

// Usuwanie zasobu z kontekstu
void deregister_resource(void *ptr);

#endif // RESOURCE_MANAGER_H
