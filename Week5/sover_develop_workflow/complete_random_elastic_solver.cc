#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> // for std::max
#include <numeric> // for std::accumulate
#include </gagarine/temporaires/zli/fftw/include/fftw3.h>
#include <omp.h>
#include </gagarine/temporaires/zli/eigen-3.4.0/Eigen/Dense>
#include <fstream> // for std::ofstream


// Function to generate meshgrid equivalent in C++
void GenerateMeshgrid(Eigen::MatrixXd& x, Eigen::MatrixXd& y, double L, int n, int m);

// Function to generate the meshgrid in Fourier space
void GenerateFrequencyMeshgrid(Eigen::MatrixXd& q_x, Eigen::MatrixXd& q_y, double L, int n, int m);

// Function to define the piecewise function phi(q)
Eigen::MatrixXd phi(const Eigen::MatrixXd& q, double phi_0, double L, double H);

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

void GenerateFrequencyMeshgrid(Eigen::MatrixXd& Q_x, Eigen::MatrixXd& Q_y, double L, int n, int m) {
   
}

// Implementation of function to define the piecewise function phi(q)
Eigen::MatrixXd phi(const Eigen::MatrixXd& q, double phi_0, double L, double H) {
    double C = phi_0; // represents phi_0
    double q_l = 2.0 * M_PI / L;
    double q_r = 2.0 * M_PI / L;
    double q_s = 2.0 * M_PI * 25.0 / L;
    double H = 0.8;

    Eigen::MatrixXd result(q.rows(), q.cols());

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < q.rows(); ++i) {
        for (int j = 0; j < q.cols(); ++j){

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