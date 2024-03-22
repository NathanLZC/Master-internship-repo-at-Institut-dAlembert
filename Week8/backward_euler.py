import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the material properties
#first we only consider one branch
'''parameters in Bugnicourt(2017) Chapter 2.2
G_0 = 2.75  # MPa
G_inf = 0.275  # MPa
G_1 = 0.3056  # MPa
'''
G_0 = 2.75  # MPa
G_1 = 2.75  # MPa
G_inf = 1/(1/G_0 + 1/G_1)  # MPa
tau_0 = 0.01  # s
eta_1 = G_1 * tau_0  # Characteristic time


#define input parameters
##time
t0 = 0
t1 = 1
dt = (t1 - t0)/1000
##load(constant)
W = 1e2  # Total load

#domain size
R = 1  # Radius of demi-sphere
L = 2  # Domain size
S = L**2  # Domain area

# Generate a 2D coordinate space
n = 300
m = 300

def contact_solver(n, m, W, S, G_t, h_profile, tol=1e-6, iter_max=10000):
    
    # Define the kernel in the Fourier domain
    q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
    q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
    QX, QY = np.meshgrid(q_x, q_y)

    kernel_fourier = np.zeros_like(QX)
    kernel_fourier = 2 / (G_t * np.sqrt(QX**2 + QY**2))
    kernel_fourier[0, 0] = 0  # Avoid division by zero at the zero frequency

    # Initial pressure distribution
    P = np.full((n, m), W / S)  # Initial guess for the pressure

    # Initialize variables for the iteration
    k = 0  # Iteration counter
    error = np.inf  # Initialize error
    h_rms = np.std(h_profile)

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

    
    while np.abs(error) > tol and k < iter_max:
        # Calculate the gradient G in the Fourier domain and transform it back to the spatial domain
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
        error = np.vdot(P, (G - np.min(G))) / (h_profile.size*h_rms*W) #/ np.linalg.norm(h_matrix)
        print(error, k, np.mean(P))
        
        k += 1  # Increment the iteration counter

    # Ensure a positive gap by updating G
    G = G - np.min(G)


    displacement_fourier = P_fourier * kernel_fourier
    displacement = np.fft.ifft2(displacement_fourier, norm='ortho').real

    return displacement, P

#######################################
###if we let k=1, we can compare the real contact area with hertz solution at t=0 and t>>\tau_0
#######################################
k = 1

alpha = G_inf + (G_1 + eta_1/dt)/(1 + G_1/G_0 + eta_1/G_0/dt)
#print(alpha)

'''
G_1_i = G_1
G_0_i = G_0
eta_i = eta_1
for i in range(0, k):
    alpha += (G_1_i + eta_i/dt)/(1 + G_1_i/G_0_i + eta_i/G_0_i/dt)
alpha += G_inf
print(alpha)
'''

beta = (eta_1/dt)/(1+G_1/G_0+eta_1/G_0/dt)
#print(beta)

'''
G_1_i = G_1
G_0_i = G_0
eta_i = eta_1
for i in range(0, k):
    beta += (eta_i/dt)/(1+G_1_i/G_0_i+eta_i/G_0_i/dt)
print(beta)
'''

gamma = (eta_1/G_0/dt)/(1+G_1/G_0+eta_1/G_0/dt)
#print(gamma)

'''
G_1_i = G_1
G_0_i = G_0
eta_i = eta_1
gamma = (eta_i/G_0_i/dt)/(1+G_1_i/G_0_i+eta_i/G_0_i/dt)
print(gamma)
'''

Surface = np.loadtxt("surface.dat")

U = np.zeros((n, m))
M = np.zeros((n, m))

for t in range(t0, t1, dt):
    #main step0: Update the effective modulus
    #effictive modulus
    G_t = G_inf + G_1 * np.exp(-t/tau_0)

    #E_star = G_t
    E_star = 1.0#clarify to be constant on both sides to update the surface profile

    #main step1: Update the surface profile
    H_new = alpha*Surface + beta*U + gamma*M

    #main step2: Update the displacement field
    U, P = contact_solver(n, m, W, S, G_t, H_new, tol=1e-6, iter_max=10000)
    ###????? how to update the loading field W_new


    A_c = np.mean(P > 0)

    #main step3: Update the partial displacement field
    M = (G_0*dt/((G_0+G_1)*dt+eta_1))*(eta_1*M /G_0/dt + (G_1+eta_1/dt)*U -eta_1*U/dt)

displacement = U 


