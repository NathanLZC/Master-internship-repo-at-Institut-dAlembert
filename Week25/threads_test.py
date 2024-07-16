import tamaas as tm
import tamaas.utils
import matplotlib.pyplot as plt
import numpy as np
import time

# parallel computation
tm.initialize(8)

# Initialize threads and fftw
tm.set_log_level(tm.LogLevel.info)  # Show progression of solver

# discretization
n = 256

# Generating surface
surface = np.load("surface.npy")
plt.imshow(surface)
plt.title("Surface Plot")
plt.colorbar()
plt.show()


# Domain size
L = 1.

model = tm.Model(tm.model_type.basic_2d, [L, L], [n, n])

model.E = 3
model.nu = 0.5
model.E_star

shear_modulus = [3] * 50

characteristic_time = np.logspace(-3, 5, 50)


# time
t0 = 0
t1 = 50
time_steps = 300
dt = (t1 - t0) / time_steps

# Const loading
W = 5e0

#solver
solver = tm.MaxwellViscoelastic(model, surface, 1e-12, dt, shear_modulus, characteristic_time)

Ac_tamaas = []

#solve for target pressure
p_target = W / (L**2) #avarge pressure
#solver.solve(p_target)

#reset the solver to avoid history accumulation
solver.reset()


start_time = time.perf_counter()

for t in np.linspace(t0, t1, time_steps):
    solver.solve(p_target)
    Ac_tamaas.append(tm.Statistics2D.contact(model.traction))

end_time = time.perf_counter()
execution_time = end_time - start_time
print("Execution time:", execution_time)
np.save("Ac_tamaas.npy", Ac_tamaas)


##### Curve fitting to get \beta

from scipy.optimize import curve_fit


Ac_tamaas = np.load("Ac_tamaas.npy")
origin_contact_area = np.load("Ac_tamaas.npy")
# Define the target function
def target_func(t, beta, t_star):
    return beta * (t - t_star)/np.min(characteristic_time)

# Extract the relevant data
x_data_1 = np.arange(t0, t1, dt)[:150] + dt
x_data_1 = np.log(x_data_1)
y_data_1 = origin_contact_area[:150]

# Perform the curve fitting
popt_1, pcov_1 = curve_fit(target_func, x_data_1, y_data_1, p0=(1,np.min(characteristic_time)))

# Extract the fitted parameters
beta_1_fit, t_star_1_fit = popt_1

# Print the fitted parameter
print("Fitted beta:", beta_1_fit)
print(pcov_1)
print("contact area percentage", np.max(Ac_tamaas)/L**2)

plt.plot(np.arange(t0, t1, dt)+dt, Ac_tamaas)
plt.plot(np.arange(t0, t1, dt)+dt, target_func(np.log(np.arange(t0, t1, dt)+dt), beta_1_fit, t_star_1_fit))
plt.xlabel("Time(s)")
plt.ylabel("Contact area($m^2$)")
plt.xscale("log")
plt.legend()
plt.title("Contact Area Over Time")
plt.show()
