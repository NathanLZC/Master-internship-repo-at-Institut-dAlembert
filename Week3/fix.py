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

p_bar = W / S  # initial load, actuall


### i am not sure if this is the right way to define the fourier 
def pressure_discritization(x, y, x0, y0):
    return p_bar
p_initial = pressure_discritization(x, y, x0, y0)# initial pressure field
p0_fourier = np.fft.fft2(p_initial, norm='ortho')  # Fourier transform of the uniform pressure

#the following are critial

#we define the frequency with q_x and q_y
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)

kernel_fourier = np.zeros_like(QX)  # Initialize the kernel array
#non_zero_indices = (QX**2 + QY**2) != 0  # Avoid division by zero
kernel_fourier = 2 / (E_star * np.sqrt(QX**2 + QY**2))
kernel_fourier[0, 0] = 0  # Set the zero frequency component to zero

# define the complementary energy as cost function
def cost_function(displacement_real, h, P):
    return np.sum(displacement_real*P)/2 - np.sum(h*P)

# define the gradient of the cost function as jacobian
def gradient(displacement_real, h):
    return displacement_real-h

jac = gradient(displacement_real, h)

 # define the length of the pressure vector
n = len(p_bar)  # p_bar is initial guess for the pressure
# first non_negative_constraint for the pressure
non_negative_constraint = minimize.LinearConstraint(np.eye(n), np.zeros(n), np.inf*np.ones(n))

# second avarage_constraint for the pressure
average_pressure_constraint = minimize.LinearConstraint(np.ones((1, n))/S, [p_bar*S], [p_bar*S])


result = minimize(cost_function, p_bar, method='trust-constr', jac=jac, constraints=[non_negative_constraint, average_pressure_constraint])




###try this

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

# Assuming you have kernel_fourier and h_matrix prepared
# Flatten kernel_fourier and h_matrix if they are not already flattened, depending on how you pass them to the function
kernel_fourier_flattened = kernel_fourier.flatten()
h_matrix_flattened = h_matrix.flatten()

# Define the initial guess for P in the Fourier domain and flatten it for optimization
p_initial_fourier_flattened = np.fft.fft2(p_initial.reshape((n, m)), norm='ortho').flatten()

# Optimization call
result = minimize(cost_function, p_initial_fourier_flattened, args=(kernel_fourier_flattened, h_matrix_flattened), method='trust-constr', options={'disp': True})

# Process the result
P_optimized_fourier = result.x.reshape((n, m))
P_optimized = np.fft.ifft2(P_optimized_fourier, norm='ortho').real

