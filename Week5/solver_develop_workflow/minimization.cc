#include <fftw3.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> // for std::max
#include <numeric> // for std::accumulate
#include <omp.h>
#include <Eigen/Dense>

// Define constants
const double W = 1e2; // Total load
const double R = 1.0; // Radius of demi-sphere
const double L = 2.0; // Domain size
const double S = L * L; // Domain area
const double E = 1e3; // Young's modulus
const double nu = 0.3; // Poisson's ratio
const double E_star = E / (1 - nu * nu); // Plane strain modulus
const int n = 100;
const int m = 100;
const double tol = 1e-6; // Tolerance for convergence
const int iter_max = 10000; // Maximum number of iterations

// Function prototypes
void generateMesh(double* x, double* y, int n, int m, double L);
void updatePressure(fftw_complex* P, fftw_complex* G, double* h_matrix, int n, int m, double alpha_0);
double alphavalue(double* P, double* G, double alpha);
template<typename T>
double findAlpha0(std::vector<T>& P, double W, double alpha_l, double alpha_r, double tol);

int main() {
    // Allocate memory for mesh grids, pressure, and displacement arrays
    double* x = new double[n * m];
    double* y = new double[n * m];
    double* h_matrix = new double[n * m];
    fftw_complex* P = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * m);
    fftw_complex* G = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * m);

    // Generate mesh grid
    generateMesh(x, y, n, m, L);

    // Define the separation h_matrix and initial pressure distribution P
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int index = i * m + j;
            double x0 = 1.0, y0 = 1.0;
            h_matrix[index] = -((x[index] - x0) * (x[index] - x0) + (y[index] - y0) * (y[index] - y0)) / (2 * R);
            P[index][0] = W / S; // Real part
            P[index][1] = 0.0; // Imaginary part
        }
    }

    //Generate the kernel



    // FFTW plans
    fftw_plan p_forward = fftw_plan_dft_2d(n, m, P, P, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan p_backward = fftw_plan_dft_2d(n, m, G, G, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Main iteration loop
    int k = 0;
    double error = std::numeric_limits<double>::infinity();
    while (std::abs(error) > tol && k < iter_max) {
        // FFT of P
        fftw_execute(p_forward);

        // Apply the kernel in Fourier domain and perform inverse FFT
        // Note: Actual implementation of multiplying with the kernel and handling of G is omitted for brevity
        fftw_execute(p_backward);



        // Update P and ensure non-negativity
        // Note: The function updatePressure should adjust P based on G and other criteria
        double alpha_0 = findAlpha0(P , W, S, -1e2, 1e2, tol, n, m);
        updatePressure(P, G, h_matrix, n, m, alpha_0);

        // Increment the iteration counter
        ++k;
    }

    // Cleanup
    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);
    fftw_free(P);
    fftw_free(G);
    delete[] x;
    delete[] y;
    delete[] h_matrix;

    // Further code for displaying results or output is omitted

    return 0;
}

void generateMesh(double* x, double* y, int n, int m, double L) {
    // Implementation of mesh grid generation is omitted for brevity

    double dx = L / (n - 1);
    double dy = L / (m - 1);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int index = i * m + j;
            x[index] = i * dx;
            y[index] = j * dy;
        }
    }

}

void updatePressure(fftw_complex* P, fftw_complex* G, double alpha_0) {
    // Implementation of pressure update is omitted for brevity


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

template<typename T>
double alphavalue(std::vector<T>& P, double W, double alpha) {
    // Implementation of finding alpha value is omitted for brevity     
    return mean(P, alpha) - W;
}

template<typename T>
double findAlpha0(std::vector<T>& P, double W, double alpha_l, double alpha_r, double tol) {

    // Expanding the search range if alpha_l and alpha_r do not bound a root
    while (sign(alphavalue(alpha_l)) == sign(alphavalue(alpha_r))) {
        alpha_r *= 2;
        // Optionally, you could throw an exception or handle the error if bounds are incorrect
    }

    // Midpoint
    double alpha_c = (alpha_l + alpha_r) / 2.0;

    // Checking if the midpoint satisfies the tolerance condition
    if (std::abs(alphavalue(alpha_c)) < tol) {
        return alpha_c;
    } else if (sign(alphavalue(alpha_l)) == sign(alphavalue(alpha_c))) {
        return findAlpha0(P, W, alpha_c, alpha_r, tol); // Narrowing the search to the right half
    } else {
        return findAlpha0(P, W, alpha_l, alpha_c, tol); // Narrowing the search to the left half
    }


    return 0.0;
}
