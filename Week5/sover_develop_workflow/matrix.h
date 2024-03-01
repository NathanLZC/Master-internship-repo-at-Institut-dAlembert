//Guard block
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

template<typename T>
class matrix {
    private:
        T* data;
        int rows;
        int cols;
    public:
        matrix(int rows, int cols);
        ~matrix();
        T& operator()(int i, int j);
        T operator()(int i, int j) const;
        int getRows() const;
        int getCols() const;
};

//implementation of the class
template<typename T>
matrix<T>::matrix(int rows, int cols) : rows(rows_), cols(cols_) {
    data = new T[rows * cols];
}

template<typename T>
matrix<T>::~matrix() {
    delete[] data;
}

template<typename T>
T& matrix<T>::operator()(int i, int j) {
    int index = i * cols + j;
    return data[index];
}

template<typename T>
T matrix<T>::operator()(int i, int j) const {
    int index = i * cols + j;
    return data[index];
}

template<typename T>
int matrix<T>::getRows() const {
    return rows;
}

template<typename T>
int matrix<T>::getCols() const {
    return cols;
}

#endif // MATRIX_H