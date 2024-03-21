import numpy as np
import matplotlib.pyplot as plt


# Define the parameters
W = 1e2  # Total load
#R = 1  # Radius of demi-sphere
L = 2  # Domain size
S = L**2  # Domain area

# Material parameters
E = 1e3  # Young's modulus
nu = 0.3  # Poisson's ratio
E_star = E / (1 - nu**2)  # Plane strain modulus

# Generate a 2D coordinate space
n = 300
m = 300
#x, y = np.meshgrid(np.linspace(0, L, n, endpoint=False), np.linspace(0, L, m, endpoint=False))

x0 = 1
y0 = 1

#################################################
#####Generate random rough surface ##############
#################################################

# Constants for the piecewise function
C = 1 # represents phi_0
q_l = 2*np.pi/L
q_r = 2*np.pi/L
q_s = 2*np.pi*25/L
H = 0.8 # Hurst exponent, represents the roughness of the surface

# Defining the piecewise function
def phi(q):
    # Initialize the result with zeros
    result = np.zeros_like(q)
    
    # Apply conditions
    mask1 = (q_l <= q) & (q < q_r)  # Condition for C
    mask2 = (q_r <= q) & (q < q_s)  # Condition for the power-law decay
    
    result[mask1] = C
    result[mask2] = C * (q[mask2] / q_r) ** (-2 * (H + 1))
    
    return result

#we define the frequency with q_x and q_y
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)

q_values = np.sqrt(QX**2 + QY**2)

phi_values = phi(q_values)

# Generate random phase
gen = np.random.default_rng()
#theta = gen.uniform(0, 2*np.pi, size=phi_values.shape)

# Generate white noise and apply PSD and phase
white_noise = gen.normal(size=phi_values.shape)
fft_noise = np.fft.fft2(white_noise)
filtered_noise = fft_noise * np.sqrt(phi_values) #* np.exp(1j * theta)
surface = np.fft.ifft2(filtered_noise).real*np.sqrt(n*m)

# Plotting
plt.figure(figsize=(10, 8))
plt.imshow(surface, cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar()
plt.title("Synthesized Rough Surface")
#plt.show()

#######################################################################
###Compute pressure and displacement with generated surface ###########
#######################################################################

h_profile = surface

# Initial pressure distribution
P = np.full((n, m), W / S)  # Initial guess for the pressure

# Define the kernel in the Fourier domain
kernel_fourier = np.zeros_like(QX)
kernel_fourier = 2 / (E_star * np.sqrt(QX**2 + QY**2))
kernel_fourier[0, 0] = 0  # Avoid division by zero at the zero frequency

# Initialize variables for the iteration
tol = 1e-6  # Tolerance for convergence
iter_max = 10000  # Maximum number of iterations
k = 0  # Iteration counter
error = np.inf  # Initialize error


def find_alpha_0(P, W, alpha_l, alpha_r, tol):
    # approximates a root, R, of f bounded 
    # by a and b to within tolerance 
    # | f(m) | < tol with m the midpoint 
    # between a and b Recursive implementation
    
    # check if a and b bound a root
    def f(alpha):
        #return np.mean(P + alpha) - W/S
        return np.mean(np.maximum(P + alpha, 0)) - W/S

    while np.sign(f(alpha_l)) == np.sign(f(alpha_r)):
        alpha_r *= 2
        #raise Exception(f"The scalars alpha_l and alpha_r do not bound a root {np.sign(f(alpha_l))} {np.sign(f(alpha_r))}")
    # get midpoint
    alpha_c = (alpha_l + alpha_r)/2
    

    if np.abs(f(alpha_c)) < tol:
        # stopping condition, report alpha_c as root
        return alpha_c
    elif np.sign(f(alpha_l)) == np.sign(f(alpha_c)):
        # case where m is an improvement on a. 
        # Make recursive call with a = m
        return find_alpha_0(P, W, alpha_c, alpha_r, tol)
    elif np.sign(f(alpha_r)) == np.sign(f(alpha_c)):
        # case where m is an improvement on b. 
        # Make recursive call with b = m
        return find_alpha_0(P, W, alpha_l, alpha_c, tol)

#RMS height(Root Mean Square height) of the surface, often denoted as R_q or S_q
h_rms = np.std(h_profile)

while np.abs(error) > tol and k < iter_max:
    # Calculate the gap G(as Gradient, see Lucas(2020)) in the Fourier domain and transform it back to the spatial domain
    P_fourier = np.fft.fft2(P, norm='ortho')
    G_fourier = P_fourier * kernel_fourier
    G = np.fft.ifft2(G_fourier, norm='ortho').real - h_profile


    # Update P by subtracting G
    P = P - G
    
    # Ensure P is non-negative
    #P = np.maximum(P, 0)
    
    # Adjust P to satisfy the total load constraint
    alpha_0 = find_alpha_0(P, W, -np.max(P), W, tol)
    #alpha_0 = find_alpha_0(P, W, -1e2, 1e2, 1e-6)
    P += alpha_0                                               # We update the pressure field inside find_alpha_0 function
    P[P < 0] = 0
    
    # Calculate the error for convergence checking
    error = np.vdot(P, (G - np.min(G))) / (surface.size*h_rms*W) #/ np.linalg.norm(h_matrix)
    print(error, k, np.mean(P))
    
    k += 1  # Increment the iteration counter


# Ensure a positive gap by updating G
G = G - np.min(G)


displacement_fourier = P_fourier * kernel_fourier
displacement = np.fft.ifft2(displacement_fourier, norm='ortho').real

# Plot the pressure field
plt.figure(figsize=(10, 8))
plt.imshow(P, cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Pressure (P)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Pressure Field')

# Plot the displacement field
plt.figure(figsize=(10, 8))
plt.imshow(displacement, cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Displacement (u_z)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement Field')
plt.show()