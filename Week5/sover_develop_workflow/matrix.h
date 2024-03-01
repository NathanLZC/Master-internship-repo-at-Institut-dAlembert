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







#endif // MATRIX_H