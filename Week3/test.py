import numpy as np
from scipy.optimize import minimize

# Define the parameters
F = 1e2  # Total force
R = 1  # Radius of demi-sphere
L = 2  # Domain size
S = L**2  # Domain area
p_bar = F / S  # Average pressure

# Define the gap function
def h(x, y, R=1):
    r = np.sqrt(x**2 + y**2)
    return -r**2 / (2 * R)

# Define the complementary energy function
def E_c(p, x, y, R=1):
    # Assuming u_z is directly related to p for simplification
    # In a real scenario, u_z should be derived from the displacement field
    u_z = p  # This is a simplification
    gap = u_z - h(x, y, R)
    return 0.5 * np.sum(u_z * p) - np.sum(p * gap)

# Objective function for minimization
def objective(p, x, y):
    return E_c(p, x, y)

# Constraint: average pressure
def constraint_avg_pressure(p):
    return np.sum(p) / S - p_bar

# Constraints
cons = [{'type': 'eq', 'fun': constraint_avg_pressure}]


# Discretize the domain
x = np.linspace(-L/2, L/2, int(np.sqrt(S)))
y = np.linspace(-L/2, L/2, int(np.sqrt(S)))
xx, yy = np.meshgrid(x, y)
xx = xx.flatten()
yy = yy.flatten()

# Assuming x and y are 1D arrays representing the discretized domain
# Ensure p_initial matches the number of points in the domain
p_initial = np.full((len(xx),), p_bar / S)  # Adjusted to match the flattened xx and yy

# Adjust E_c to correctly handle shapes
def E_c(p, x, y, R=1):
    u_z = p  # Direct relationship for simplification

    # Calculate the gap for each point in the domain
    gap = u_z - h(x, y, R)  # Ensure h(x, y, R) returns an array matching u_z in shape

    # Calculate complementary energy
    return 0.5 * np.sum(u_z * p) - np.sum(p * gap)

# Perform the minimization
result = minimize(lambda p: objective(p, xx, yy), p_initial, method='SLSQP', constraints=cons, options={'disp': True})

if result.success:
    optimized_p = result.x
    print("Optimization successful.")
else:
    print("Optimization failed.")
