import numpy as np
from scipy.optimize import minimize

#define the parameters
W = 1e2  # total load
R = 1  # radius of demi-sphere
L = 2  # domain size
S = L**2  # domain area



#define material parameters
E = 1e3  # Young's modulus
nu = 0.3  # Poisson's ratio
E_star = E / (1 - nu**2)  # plane strain modulus 

#We generate a 2D coordinate space
n = 100
m = 100

x, y = np.meshgrid(np.linspace(0, L, n, endpoint=False), np.linspace(0, L, m, endpoint=False))#notice here that we use endpoint=False to avoid having the last point

x0 = 1
y0 = 1

# Define the separation h          #Since we are using the same model, we can use the same h function
h = - (x**2 + y**2)/(2*R)

p_bar = W / S  

h_norm = np.abs(h(x, y, R)) #norm of the h function

alpha_0 = 1e-2  # learning rate
tol = 1e-6  # tolerance for convergence
iter_max = 1000  # maximum number of iterations

# Initialize P with the average constant load initial guess
P = np.full((n, m), p_bar)


P = p_bar #Initial guess for the pressure

for k in range(0, iter_max): 
    if np.abs(error) < tol:

        #gradient(gap)
        G = 
        P = P - G




        p = 

        error =



        k += 1
        if np.abs(error) < tol:
            p_bar = p
        else:
            print('Convergence failed')

    G = G - min(0, G)  # project the gradient onto the non-negative orthant
            



