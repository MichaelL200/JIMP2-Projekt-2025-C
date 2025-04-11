#ifndef EIGENVECTORS_H
#define EIGENVECTORS_H

// Struktura do przechowywania zmiennych potrzebnych do metody Lanczosa oraz do obliczenia przybliżeń wektorów własnych macierzy Laplace'a grafu
typedef struct
{
    // rozmiar macierzy Laplace'a
    int n;
    // liczba iteracji
    int m;
    // baza ortonormalna Q (przybliżenia wektorów własnych)
    double* Q;
    // macierz trójdiagonalna T
    double* alpha;
    double* beta;
    // wartości własne macierzy T
    double *eigenvalues_T; // theta
    // wektory własne macierzy T (macierz)
    double *eigenvectors_T; // y
    // przybliżone wektory własne macierzy L (macierz)
    double *eigenvectors_L; // x
} LanczosEigenV;

// Funkcja obliczająca iloczyn skalarny dwóch wektorów o długości n
double dot_product(double *v1, double *v2, int n);

// Funkcja obliczająca normę euklidesową wektora o długości n
double norm(double *v, int n);

// Funkcja mnożąca macierz M (n x n, przechowywana w porządku wierszowym) przez wektor v
// Wynik zapisywany jest w tablicy result
void mat_vec_multiply(double* M, double* v, double* result, int n);

// Funkcja testująca poprzednie
void test_ev();

#endif // EIGENVECTORS_H