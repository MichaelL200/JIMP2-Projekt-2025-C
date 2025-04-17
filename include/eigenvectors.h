#ifndef EIGENVECTORS_H
#define EIGENVECTORS_H

// Struktura do przechowywania zmiennych potrzebnych do metody Lanczosa oraz do obliczenia przybliżeń wektorów własnych macierzy Laplace'a grafu
typedef struct
{
    // rozmiar macierzy Laplace'a
    int n;
    // liczba iteracji
    int m;
    // baza ortonormalna V z wektorami v (baza przestrzeni Kryłowa)
    double* V;
    // wektory w (macierz)
    double* W;
    // macierz trójdiagonalna T
    double* alpha;
    double* beta;
    // wartości własne macierzy T
    double* theta;
    // wektory własne macierzy T (macierz)
    double* Y;
    // przybliżone wektory własne macierzy L (macierz)
    double* X;
} LanczosEigenV;
// Nazwy takie same jak na Wikipedii - Lanczos algorithm
// https://en.wikipedia.org/wiki/Lanczos_algorithm

// Funkcja obliczająca iloczyn skalarny dwóch wektorów o długości n
double dot_product(double *v1, double *v2, int n);

// Funkcja obliczająca normę euklidesową wektora o długości n
double norm(double *v, int n);

// Funkcja mnożąca macierz M (n x n, przechowywana w porządku wierszowym) przez wektor v. Wynik zapisywany w tablicy result
void mat_vec_multiply_d(double* M, double* v, double* result, int n);
void mat_vec_multiply_i(int* M, double* v, double* result, int n);

// Funkcja testująca funkcje: dot_product, norm i mat_vec_multiply
void test1();

// Funkcja ortogonalizująca wektor w stosunku do poprzednich wektorów
void orthogonalize(double *v, double *V, int n, int j);

// Funkcja normalizująca wektor
int normalize(double *v, int n);

// Funkcja testująca funkcje orthogonalize i normalize
void test2();

// Inicjalizacja wartości obiektu struktury
void lanczos_init(LanczosEigenV *l, int n, int m);

// Losowanie dowolnego wektora v₁
void lanczos_v1_init(LanczosEigenV *l);

// Pierwszy inicjalizacyjny krok iteracyjny metody Lanczosa
void lanczos_initial_step(LanczosEigenV *l, int* A);

// Iteracje metody Lanczosa dla j = 2, ..., m
void lanczos(LanczosEigenV *l, int *A);

// Funkcja testująca funkcje metody Lanczosa
void test3();

#endif // EIGENVECTORS_H
