# Midterm Script Part 2
import numpy as np
import matplotlib.pyplot as plt

'''Constants'''
L = 0.04 # m; length of fin
t = 0.001 # m; thickness of fin
M = 16 # number of fins 
b = 0.006 # m; base plate's thickness 
w = 0.06 # m; base plate's length and width 
g = 0.00275 # m; air gap between each fin
q = 250 # W; electrical component power generation 
h_list = np.arange(10, 110, 10) # W/m^2K; convective coefficent
n = 10 # number of nodes the fin is broken down to
dx = L / n # m
dt = 0.01 # s; time step
T_inf = 25 # oC; ambient temperature

# Aluminum
k_al = 200 # W/mK
rho_al = 2700 # kg/m^3
c_al = 910 # J/kgK

'''Setup'''
T_list = [] # List of temperature arrays

'''Simulation'''
it = 0 # iterator
for h in h_list:
    k = k_al
    rho = rho_al
    c = c_al
    
    Bi = ((h * t / 4) / k)
    R = np.array([1 / (h * w * g / 2), # 0, inf
                  (dx / 2) / (k * w * t / 2), # 0, 1
                  dx / (k * w * t / 2), # i, i + 1
                  (1 / (h * dx)) * (1 + Bi), # i, inf
                  1 / ((1 / ((1 / (h * dx)) + ((t / 4) / (k * dx)))) + (1 / ((1 / (h * t / 2)) + ((dx / 2) / (k * t / 2)))))]) # n, inf
    C = np.ones(n + 1) * rho * c * w * (t / 2) * dx # J/oC; capacitance
    C[0] = rho * c * w * b * (0.5 * (g + t)) # J/oC; capacitance of baseplate
    
    tau = np.ones(n + 1) * C[2] / ((1 / R[2]) + (1 / R[2]) + (1 / (R[3] * (1 + Bi))))
    tau[1] = C[1] / ((1 / R[1]) + (1 / R[2]) + (1 / (R[3] * (1 + Bi))))
    tau[0] = C[0] / ((1 / R[1]) + (1 / R[0]))
    tau[n] = C[n] / ((1 / R[2]) + (1 / (R[4] * (1 + Bi))))
    
    time =  np.arange(0, 500, dt) # s
    T = np.zeros((n + 1, len(time))) # oC
                 
    p = 0
    
    # Power on component until steady state is reached
    q = 250 # W; electrical component power generation
    while (time[p] < time[-1] / 2):
        # Base Plate
        T[0, p + 1] = T[0, p] + (dt / C[0]) * ((q / (2 * M)) + ((T[1, p] - T[0, p]) / R[1]) + ((T_inf - T[0, p]) / R[0]))
        
        # Node i = 1 to i = n - 1
        for i in range(1, n):
            if (i == 1):
                T[i, p + 1] = T[i, p] + (dt / C[i]) * (((T[i - 1, p] - T[i, p]) / R[1]) + ((T[i + 1, p] - T[i, p]) / R[2]) + ((T_inf - T[i, p]) / (R[3] * (1 + Bi)))) # oC
            else:
                T[i, p + 1] = T[i, p] + (dt / C[i]) * (((T[i - 1, p] - T[i, p]) / R[2]) + ((T[i + 1, p] - T[i, p]) / R[2]) + ((T_inf - T[i, p]) / (R[3] * (1 + Bi)))) # oC
        
        # Node n
        T[n, p + 1] = T[n, p] + (dt / C[n]) * (((T[n - 1, p] - T[n, p]) / R[2]) + ((T_inf - T[n, p]) / ((R[4]) * (1 + Bi)))) # oC
    
        p = p + 1
    
    # Turn off component and cool until ambient temperature is reached    
    q = 0
    while (time[p] < time[-1]):
        # Base Plate
        T[0, p + 1] = T[0, p] + (dt / C[0]) * ((q / (2 * M)) + ((T[1, p] - T[0, p]) / R[1]) + ((T_inf - T[0, p]) / R[0]))
        
        # Node i = 1 to i = n - 1
        for i in range(1, n):
            if (i == 1):
                T[i, p + 1] = T[i, p] + (dt / C[i]) * (((T[i - 1, p] - T[i, p]) / R[1]) + ((T[i + 1, p] - T[i, p]) / R[2]) + ((T_inf - T[i, p]) / (R[3] * (1 + Bi)))) # oC
            else:
                T[i, p + 1] = T[i, p] + (dt / C[i]) * (((T[i - 1, p] - T[i, p]) / R[2]) + ((T[i + 1, p] - T[i, p]) / R[2]) + ((T_inf - T[i, p]) / (R[3] * (1 + Bi)))) # oC
        
        # Node n
        T[n, p + 1] = T[n, p] + (dt / C[n]) * (((T[n - 1, p] - T[n, p]) / R[2]) + ((T_inf - T[n, p]) / ((R[4]) * (1 + Bi)))) # oC
    
        p = p + 1
        
    T_list.append(T)
    
'''Plots'''
# Base Temperature vs. Convective Coefficient at Steady State
T_list_inc = []
for Temperature in T_list:
    T_list_inc.append(Temperature[0, 25000])
plt.figure()
plt.plot(h_list, T_list_inc)
plt.title('Base Temperature vs. Convective Coefficient (at Steady-State) for Aluminum')
plt.xlabel('Convective Coefficient [W/m^2/K]')
plt.ylabel('Base Temperature [oC]')
plt.show()

# Temperature vs. Position for various convective coefficients
x_0 = np.array([b / 2])
x = np.arange(b + (dx / 2), (L + b) + (dx / 2), dx) # m
x = np.concatenate((x_0, x), axis = 0)
plt.figure()
i = 0
for Temperature in T_list:
    plt.plot(x, Temperature[:, 25000], label = 'h =' + str(h_list[i]) + ' W/m^2/K')
    i = i + 1
plt.title('Temperature vs. Position (at Steady-State) for Various Convective Coefficients (Aluminum)')
plt.xlabel('Position [m]')
plt.ylabel('Temperature [oC]')
plt.legend()
plt.show()

