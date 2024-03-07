//guard block
#ifndef ELASTIC_H
#define ELASTIC_H

#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <omp.h>

template<typename T>
void generateMesh(T* x, T* y, T* h_matrix, int n, int m, T L);

template<typename T>
void generateMesh(T* x, T* y, T* h_matrix, int n, int m, T L) {
    T dx = L / (n - 1);
    T dy = L / (m - 1);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int index = i * m + j;
            x[index] = i * dx;
            y[index] = j * dy;
            h_matrix[index] = 0.0;
        }
    }
}

//#include "mesh.tpp"//why this template is not working

#endif // ELASTIC_H