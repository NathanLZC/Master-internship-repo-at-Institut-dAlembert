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

# Partial displacement
M = 1e-3  # m


#######################################
###Q1: if we let k=1, what result do we espect?
#######################################
###Q2: for a general case, I should define alpha, beta, gamma function and use
###    discretization points as the value of k.
#######################################

alpha = G_inf + (G_1 + eta_1/dt)/(1 + G_1/G_0 + eta_1/G_0/dt)

beta = (eta_1/dt)/(1+G_1/G_0+eta_1/G_0/dt)

gamma = (eta_1/G_0/dt)/(1+G_1/G_0+eta_1/G_0/dt)


H_old = 1

def solve_U()

for t in np.arange(t0, t1, dt):

    H_new = alpha*H_old + beta*W + gamma*M


'''
#surface profile
h_profile =  np.loadtxt("surface.dat")
'''







