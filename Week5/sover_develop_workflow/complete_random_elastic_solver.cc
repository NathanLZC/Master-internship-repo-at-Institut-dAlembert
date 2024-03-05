#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> // for std::max
#include <numeric> // for std::accumulate
#include </gagarine/temporaires/zli/fftw/include/fftw3.h>
#include <omp.h>
#include </gagarine/temporaires/zli/eigen-3.4.0/Eigen/Dense>
#include <random> // for std::random_device, std::mt19937, std::normal_distribution
#include <fstream> // for std::ofstream


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
double mean(std::vector<T>& v, double n);

// Function to compute the sign of a vector
template<typename T>
std::vector<int> sign(const std::vector<T>& v);

// Function to find the alpha value for steepest descent algorithm energy minimization
template<typename T>
double findAlpha0(std::vector<T>& P, double W, double alpha_l, double alpha_r, double tol);

// Generate kernel function
Eigen::MatrixXd generateKernel(const Eigen::MatrixXd& QX, const Eigen::MatrixXd& QY, double E_star);

// Function to minimize the energy
Eigen::MatrixXd computeDisplacment(const Eigen::MatrixXd& surface, Eigen::MatrixXd& P, Eigen::MatrixXd kernel_fourier, double W, double S, double error, double tol, int k, int maxIter);


////////////////////////////////////////
//implementation of the functions
////////////////////////////////////////

void GenerateMeshgrid(Eigen::VectorXd& x, Eigen::VectorXd& y,Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double L, int n, int m) {
    x = Eigen::VectorXd::LinSpaced(n, 0, L);
    y = Eigen::VectorXd::LinSpaced(m, 0, L);
    X = Eigen::MatrixXd::Zero(n, m);
    Y = Eigen::MatrixXd::Zero(n, m);
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            X(i, j) = x(i);
            Y(i, j) = y(j);
        }
    }
}

// Implementation of function to implement numpy.fft.fftfreq in  C++
Eigen::VectorXd fftfreq(int n, double d=1.0){
    Eigen::VectorXd result(n);
    if (n % 2 == 0) {
        int N = n / 2;
        #pragma omp parallel for
        for(int i = 0; i < N; ++i) {
            result(i) = i / (N * d);
        }
        #pragma omp parallel for
        for(int i = N; i < n; ++i) {
            result(i) = (i - n) / (N * d);
        }
    } else {
        int N = (n - 1) / 2;
        #pragma omp parallel for
        for(int i = 0; i < N; ++i) {
            result(i) = i / (N * d);
        }
        #pragma omp parallel for
        for(int i = N; i < n; ++i) {
            result(i) = (i - n) / (N * d);
        }
    }
    return result;
}

void GenerateFrequencyMeshgrid(Eigen::VectorXd& q_x, Eigen::VectorXd& q_y,Eigen::MatrixXd& Q_x, Eigen::MatrixXd& Q_y, double L, int n, int m){
    q_x = fftfreq(n, L / n);
    q_y = fftfreq(m, L / m);
    Q_x = Eigen::MatrixXd::Zero(n, m);
    Q_y = Eigen::MatrixXd::Zero(n, m);
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            Q_x(i, j) = q_x(i);
            Q_y(i, j) = q_y(j);
        }
    }
}

// Implementation of function to define the q value
Eigen::MatrixXd computeQValues(const Eigen::VectorXd& Q_x, const Eigen::VectorXd& Q_y){
    Eigen::MatrixXd result(Q_x.rows(), Q_x.cols());
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < Q_x.rows(); ++i) {
        for(int j = 0; j < Q_x.cols(); ++j) {
            result(i, j) = std::sqrt(Q_x(i, j) * Q_x(i, j) + Q_y(i, j) * Q_y(i, j));
        }
    }
    return result;
}

// Implementation of function to define the piecewise function phi(q)
Eigen::MatrixXd phi(const Eigen::MatrixXd& q, double phi_0, double q_l, double q_r, double q_s, double H) {
    double C = phi_0; // represents phi_0
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(q.rows(), q.cols());
    
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < q.rows(); ++i) {
        for (int j = 0; j < q.cols(); ++j){
            double val = q(i, j);
            if (q_l <= val && val < q_r){
                result(i, j) = C;
            } else if (q_r <= val && val < q_s){
                result(i, j) = C * std::pow((val / q_r), -2*(H+1));
            }
        }
    }
}

// Implementation to generate white noise//https://cplusplus.com/reference/random/mt19937/
Eigen::MatrixXd generateWhiteNoise(int rows, int cols) {
    std::random_device rd{}; //std::random_device is a mechanism provided by C++ to generate a non-deterministic random number. 
                             //It's often used as a seed for more complex random number generators. The {} initializes an instance of std::random_device.
    std::mt19937 gen{rd()};  //std::mt19937 is a pseudo-random number generator (PRNG) based on the Mersenne Twister algorithm. 
                             //It's known for producing high-quality random numbers and has a very long period (the sequence of random numbers it generates before repeating is extremely long).
    std::normal_distribution<> d{0,1}; // mean 0, standard deviation 1

    Eigen::MatrixXd whiteNoise(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            whiteNoise(i, j) = d(gen);
        }
    }
    return whiteNoise;
}

Eigen::MatrixXd generateRandomSurface(const Eigen::MatrixXd& phiValues, int n, int m){
    // Generate white noise using the predefined function
    Eigen::MatrixXd whiteNoise = generateWhiteNoise(n, m);

    // Allocate memory for the surface(return value)
    Eigen::MatrixXd surface = Eigen::MatrixXd::Zero(n, m);

    // Allocate memory for FFTW inputs/outputs
    //fftw_complex* fftInput = (fftw_complex*)fftw_malloc(n * m * sizeof(fftw_complex));
    fftw_complex* fftOutput = (fftw_complex*)fftw_malloc(n * (m/2 + 1) * sizeof(fftw_complex));

    // FFTW plan
    //FFT and IFFT Operations: These are performed using FFTW's r2c (real-to-complex) and c2r (complex-to-real) transformations, 
    //suitable for handling the real-valued nature of the white noise and the resulting surface.
    fftw_plan p_forward = fftw_plan_dft_r2c_2d(n, m, whiteNoise.data(), fftOutput, FFTW_ESTIMATE);

    // Excute the foward plan
    fftw_execute(p_forward);

    // Apply the filter
    #pragma omp parallel for collapse(2)
    for (int i=0; i < n; ++i) {
        for (int j=0; j < m/2 + 1; ++j) { ////??????
            int index = i * (m/2 + 1)+ j;
            double filter = std::sqrt(phiValues(i, j));
            fftOutput[index][0] *= filter; // Real part
            fftOutput[index][1] *= filter; // Imaginary part
        }
    }

    // FFTW plan for inverse transform
    fftw_plan p_backward = fftw_plan_dft_c2r_2d(n, m, fftOutput, surface.data(), FFTW_ESTIMATE);


    // Execute the backward plan
    fftw_execute(p_backward);

    // Deallocate memory
    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);
    //fftw_free(fftInput);
    fftw_free(fftOutput);

    return surface;
}

void SaveSurfaceToFile(const Eigen::MatrixXd& surface, const std::string& filename){
    std::ofstream file(filename);
    if (file.is_open()) {
        file << surface << '\n';
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
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

//inplementation of the kernel function
Eigen::MatrixXd generateKernel(const Eigen::MatrixXd& QX, const Eigen::MatrixXd& QY, double E_star){
    Eigen::MatrixXd kernel_fourier = Eigen::MatrixXd::Zero(QX.rows(), QX.cols());
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < QX.rows(); ++i) {
        for(int j = 0; j < QX.cols(); ++j) {
            if (!(i == 0 && j == 0)) {
                double qMagnitude = std::sqrt(QX(i, j) * QX(i, j) + QY(i, j) * QY(i, j));
                kernel_fourier(i, j) = 2 / (E_star * qMagnitude);
            }
        }
    }
    return kernel_fourier;
}

// Function to minimize the energy to compute the displacement field
Eigen::MatrixXd computeDisplacment(const Eigen::MatrixXd& surface, Eigen::MatrixXd& P, Eigen::MatrixXd kernel_fourier, double W, double S, double error, double tol, int k, int maxIter){
    int n = P.rows();
    int m = P.cols();
    
    // Allocate memory for the pressure field, gap function and displacement field
    fftw_complex* P_fourier = (fftw_complex*)fftw_malloc(n * (m/2+1)* sizeof(fftw_complex));
    fftw_complex* G_fourier = (fftw_complex*)fftw_malloc(n * (m/2+1)* sizeof(fftw_complex));
    fftw_complex* displacement_fourier = (fftw_complex*)fftw_malloc(n * (m/2+1)* sizeof(fftw_complex));

    // G
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(n, m);

    // Displacement
    Eigen::MatrixXd displacement = Eigen::MatrixXd::Zero(n, m);


    while (std::abs(error) > tol && k < maxIter) {
        // Perform FFT on P
        // FFTW plan
        fftw_plan p_forward = fftw_plan_dft_r2c_2d(n, m, P.data(), P_fourier, FFTW_ESTIMATE);

        // Execute the forward plan
        fftw_execute(p_forward);

        G_fourier = P_fourier * kernel_fourier;

        fftw_plan p_backward = fftw_plan_dft_c2r_2d(n, m, G_fourier, G, FFTW_ESTIMATE);
        
        // Apply kernel in Fourier domain and perform inverse FFT to get G
        // Note: Actual application of kernel_fourier and inverse FFT needed
        
        // Subtract h_profile from G
        G = G - surface;
        
        // Update P by subtracting G
        P = P - G;
        
        // Ensure P is non-negative
        P = P.cwiseMax(0.0);
        
        // Adjust P to satisfy the total load constraint
        double alpha_0 = findAlpha0(P, W / S, P.minCoeff(), W, tol);
        P += alpha_0;
        P = P.cwiseMax(0.0);
        
        // Calculate the error for convergence checking
        // This is a simplified version; adjust as needed
        error = (P.array() * (G.array() - G.minencapsulate Coeff())).sum() / (P.size() * W);
        
        std::cout << "Error: " << error << ", Iteration: " << k << std::endl;
        
        k++;  // Increment the iteration counter

        // Deallocate memory
        fftw_destroy_plan(p_forward);
        fftw_destroy_plan(p_backward);


    }


    //FFTW plan for displacement
    fftw_plan p_forward = fftw_plan_dft_r2c_2d(n, m, displacement.data(), displacement_fourier, FFTW_ESTIMATE);

    // Execute the forward plan
    fftw_execute(p_forward);

    displacement_fourier = P_fourier * kernel_fourier;

    //FFTW plan for inverse transform
    fftw_plan p_backward = fftw_plan_dft_c2r_2d(n, m, displacement_fourier, displacement.data(), FFTW_ESTIMATE);

    // Execute the backward plan
    fftw_execute(p_backward);


    fftw_free(P_fourier);
    fftw_free(G_fourier);

    return displacement;
}



////////////////////////////////////////


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


