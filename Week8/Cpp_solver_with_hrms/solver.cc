#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> // for std::max
#include <numeric> // for std::accumulate
#include </gagarine/temporaires/zli/fftw/include/fftw3.h>
#include <omp.h> //#include </opt/homebrew/Cellar/libomp/18.1.0/include/omp.h> // for Macbook
#include </gagarine/temporaires/zli/eigen-3.4.0/Eigen/Dense>
#include <random> // for std::random_device, std::mt19937, std::normal_distribution
#include <fstream> // for std::ofstream
#include "viscoelastic_contact.h"

///////////////////////////////////////////////////////
///here is the main function, we consider k=1, one-branch model 
///for the viscoelastic contact problem
///////////////////////////////////////////////////////



// main function
int main(){
    // Define the loading
    double W = 1e2;

    // Define the domain size
    double L = 2;
    double S = L * L;
    int n = 300;
    int m = 300;
    double dx = L / n;
    double dy = L / m;

    // Material parameters
    double E = 1e3; // Young's modulus
    double nu = 0.3;  // Poisson's ratio
    double E_star = E / (1 - nu * nu);  // Plane strain modulus

    // Generate meshgrid
    Eigen::VectorXd x, y;
    Eigen::MatrixXd X, Y;
    GenerateMeshgrid(x, y, X, Y, L, n, m);

    // Generate frequency meshgrid
    Eigen::VectorXd q_x, q_y;
    Eigen::MatrixXd Q_x, Q_y;
    GenerateFrequencyMeshgrid(q_x, q_y, Q_x, Q_y, L, n, m);

    // Define the q value
    Eigen::MatrixXd q = computeQValues(Q_x, Q_y);

    // Define the piecewise function phi(q)
    double phi_0 = 1.0;
    double q_l = 2*M_PI/L;
    double q_r = 2*M_PI/L;
    double q_s = 2*M_PI*25/L;
    double H = 0.8;
    Eigen::MatrixXd Phi = phi(q, phi_0, q_l, q_r, q_s, H);

    // Generate random phase and white noise
    Eigen::MatrixXd surface = generateRandomSurface(Phi, n, m);

    // Save the surface to a file
    SaveSurfaceToFile(surface, "surface.dat");

    // Initial guess for the pressure
    Eigen::MatrixXd P_initial = Eigen::MatrixXd::Constant(n, m, W / S);

    // Generate kernel function
    Eigen::MatrixXd kernel_fourier = generateKernel(Q_x, Q_y, E_star);

    // Update the pressure field
    double tol = 1e-6;
    int maxIter = 1000;       
    double error = std::numeric_limits<double>::infinity();
    int k = 0;

    //Compute the displacement field by energy minimization
    Eigen::MatrixXd displacement = computeDisplacment(surface, P_initial, kernel_fourier, W, S, error, tol, k, maxIter);

    //Save displacement to a file
    SaveSurfaceToFile(displacement, "displacement.dat");

    return 0;

}
