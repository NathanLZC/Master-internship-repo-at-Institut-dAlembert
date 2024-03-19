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
W = 1000

#domain size
R = 1  # Radius of demi-sphere
L = 2  # Domain size
S = L**2  # Domain area

# Generate a 2D coordinate space
n = 300
m = 300

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
    alpha = G_inf + (G_1_i + eta_i/dt)/(1 + G_1_i/G_0_i + eta_i/G_0_i/dt)
print(alpha)
'''

beta = (eta_1/dt)/(1+G_1/G_0+eta_1/G_0/dt)
#print(beta)

'''
G_1_i = G_1
G_0_i = G_0
eta_i = eta_1
for i in range(0, k):
    beta = (eta_i/dt)/(1+G_1_i/G_0_i+eta_i/G_0_i/dt)
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

H_old = np.loadtxt("surface.dat")



for t in range(t0, t1, dt):

    

    H_new = alpha*H_old + beta*W + gamma*M
'''









