#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> // for std::max
#include <numeric> // for std::accumulate
#include <fftw3.h>
#include <omp.h>
#include <Eigen/Dense>
#include <fstream> // for std::ofstream


// Function to generate meshgrid equivalent in C++
void GenerateMeshgrid(Eigen::MatrixXd& x, Eigen::MatrixXd& y, double L, int n, int m);

// Function to define the piecewise function phi(q)
Eigen::VectorXd phi(const Eigen::VectorXd& q, double L, double H);

// Function to generate random phase and white noise
void GenerateRandomSurface(Eigen::MatrixXd& surface, const Eigen::VectorXd& q_x, const Eigen::VectorXd& q_y, int n, int m);

// Function to save data to a file
void SaveSurfaceToFile(const Eigen::MatrixXd& surface, const std::string& filename);


void GenerateMeshgrid(Eigen::MatrixXd& x, Eigen::MatrixXd& y, double L, int n, int m) {
    x = Eigen::MatrixXd(n, m);
    y = Eigen::MatrixXd(n, m);
    double dx = L / (n - 1);
    double dy = L / (m - 1);

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            x(i, j) = i * dx;
            y(i, j) = j * dy;
        }
    }
}


//need to accelerate this function
template<typename T>
double mean(std::vector<T>& v, double n) {
    T sum = std::accumulate(v.begin(), v.end(), T(0)) + n * v.size();
    return sum / static_cast<double>(v.size());
}

//sign function
template<typename T>
std::vector<int> sign(const std::vector<T>& v) {
    std::vector<int> signs(v.size());
    for(size_t i = 0; i < v.size(); ++i) {
        signs[i] = (v[i] > T(0)) - (v[i] < T(0));
    }
    return signs;
}

void GenerateRandomSurface(Eigen::MatrixXd& surface, const Eigen::VectorXd& q_x, const Eigen::VectorXd& q_y, int n, int m){
    



}