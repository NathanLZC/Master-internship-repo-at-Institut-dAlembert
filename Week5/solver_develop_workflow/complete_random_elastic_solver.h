//Guard block
#pragma once

#include <iostream>
#include </gagarine/temporaires/zli/fftw/include/fftw3.h>
#include </gagarine/temporaires/zli/eigen-3.4.0/Eigen/Dense>

// Function to generate meshgrid equivalent in C++
void GenerateMeshgrid(Eigen::VectorXd& x, Eigen::VectorXd& y,Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double L, int n, int m);

// Function to realize numpy.fft.fftfreq in  C++
Eigen::VectorXd fftfreq(int n, double d = 1.0);

// Function to generate the meshgrid in Fourier space
void GenerateFrequencyMeshgrid(Eigen::VectorXd& q_x, Eigen::VectorXd& q_y,Eigen::MatrixXd& Q_x, Eigen::MatrixXd& Q_y, double L, int n, int m);

// Function to define the q value
Eigen::MatrixXd computeQValues(const Eigen::VectorXd& Q_x, const Eigen::VectorXd& Q_y);

// Function to define the piecewise function phi(q)
Eigen::MatrixXd phi(const Eigen::MatrixXd& q, double phi_0, double q_l, double q_r, double q_s, double H);

//Function to Generate white noise
Eigen::MatrixXd generateWhiteNoise(int rows, int cols);

// Function to generate random phase and white noise
Eigen::MatrixXd generateRandomSurface(const Eigen::MatrixXd& phiValues, int n, int m);

// Function to save data to a file
void SaveSurfaceToFile(const Eigen::MatrixXd& surface, const std::string& filename);

// Function to compute the mean of a vector
template<typename T>
double mean(const std::vector<T>& v);


// Function to compute the sign of a vector
/*
template<typename T>
std::vector<int> sign(const std::vector<T>& v);
*/

int sign(double value);

// Function to define the alpha value
double alphavalue(double alpha, const Eigen::MatrixXd& P, double W, double S);

// Function to find the alpha value for steepest descent algorithm energy minimization
double findAlpha0(const Eigen::MatrixXd& P, double W, double alpha_l, double alpha_r, double tol, double S);

// Generate kernel function
Eigen::MatrixXd generateKernel(const Eigen::MatrixXd& QX, const Eigen::MatrixXd& QY, double E_star);

// Function to minimize the energy
Eigen::MatrixXd computeDisplacment(const Eigen::MatrixXd& surface, Eigen::MatrixXd& P, Eigen::MatrixXd kernel_fourier, double W, double S, double error, double tol, int k, int maxIter);
