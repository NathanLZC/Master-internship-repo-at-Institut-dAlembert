import numpy as np
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint


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
h_matrix = - (x**2 + y**2)/(2*R)

# Initial pressure distribution
p_initial = np.full((n, m), W / S).flatten()  # Flatten for optimization

# Fourier transform of initial pressure distribution
P_fourier_flattened = np.fft.fft2(p_initial.reshape(n, m), norm='ortho').flatten()  # Keep flattened


#the following are critial

#we define the frequency with q_x and q_y
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)

kernel_fourier = np.zeros_like(QX)  # Initialize the kernel array
'''
#non_zero_indices = (QX**2 + QY**2) != 0  # Avoid division by zero
kernel_fourier = 2 / (E_star * np.sqrt(QX**2 + QY**2))
kernel_fourier[0, 0] = 0  # Set the zero frequency component to zero
'''
# ?????
Q_magnitude = np.sqrt(QX**2 + QY**2)
Q_magnitude[0, 0] = 1  # Avoid division by zero
kernel_fourier = 2 / (E_star * Q_magnitude)
kernel_fourier[0, 0] = 0  # Set the zero frequency component to zero


def cost_function(P_fourier_flattened, kernel_fourier, h_matrix):
    # Ensure kernel_fourier and h_matrix are correctly shaped and prepared before this function is called
    
    # Reshape P_fourier_flattened back to its 2D shape to perform operations
    P_fourier = P_fourier_flattened.reshape((n, m))
    
    # Calculate u_z in the Fourier domain
    u_z_fourier = P_fourier * kernel_fourier
    
    # Transform u_z back to the spatial domain to compute the cost function
    u_z = np.fft.ifft2(u_z_fourier, norm='ortho').real
    
    # Calculate the first term of the cost function: 0.5 * sum(u_z * P), where P is also in the spatial domain
    P = np.fft.ifft2(P_fourier, norm='ortho').real
    term1 = 0.5 * np.sum(u_z * P)
    
    # Calculate the second term of the cost function: sum(P * h)
    term2 = np.sum(P * h_matrix)
    
    # The cost function is the difference of the two terms
    cost = term1 - term2
    return cost

# define the gradient of the cost function as jacobian
def gradient(P_fourier_flattened, kernel_fourier, h_matrix):
    # Ensure kernel_fourier and h_matrix are correctly shaped and prepared before this function is called
    
    # Reshape P_fourier_flattened back to its 2D shape to perform operations
    P_fourier = P_fourier_flattened.reshape((n, m))
    
    # Calculate u_z in the Fourier domain
    u_z_fourier = P_fourier * kernel_fourier
    
    # Transform u_z back to the spatial domain to compute the cost function
    u_z = np.fft.ifft2(u_z_fourier, norm='ortho').real

    return u_z - h_matrix




# Flatten kernel_fourier and h_matrix if they are not already flattened, depending on how you pass them to the function
kernel_fourier_flattened = kernel_fourier.flatten()
h_matrix_flattened = h_matrix.flatten()


 # define the length of the pressure vector
n = len(P_fourier_flattened)  # p_bar is initial guess for the pressure
# first non_negative_constraint for the pressure
non_negative_constraint = LinearConstraint(np.eye(n), np.zeros(n), np.inf*np.ones(n))

# second avarage_constraint for the pressure
average_pressure_constraint = LinearConstraint(np.ones((1, n))/S, [P_fourier_flattened*S], [P_fourier_flattened*S])

jac = gradient(P_fourier_flattened, kernel_fourier_flattened, h_matrix)



result = minimize(cost_function, P_fourier_flattened, method='trust-constr', jac=jac, constraints=[non_negative_constraint, average_pressure_constraint])


# Process the result
# .x refers to an attribute of the object result, which is returned by the minimize function from the scipy.optimize module in Python.
P_optimized_fourier = result.x.reshape((n, m))

displacement_optimized_fourier = P_optimized_fourier * kernel_fourier
displacement_optimized_real = np.fft.ifft2(displacement_optimized_fourier, norm='ortho').real

print(displacement_optimized_real)