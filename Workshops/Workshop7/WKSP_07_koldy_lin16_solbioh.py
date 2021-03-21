import numpy as np
import matplotlib.pyplot as plt

'''USER SETUP'''
n = int(input('Please input a number of nodes: \n')) # nummber of nodes
dt = float(input('Please input a dt value: \n')) # s; timestep

'''CONSTANTS'''
# Enviornment
T_air = 30 # oC; air temperature
T_sky = 0 # oC; sunny temperature
T_cloud = 10 # oC; cloudy temperature

# Asphalt Pavement
L = 0.1 # m; length
A = 1 # m^2; area
k = 0.75 # W/mK
rho = 2360 # kg/m^3; density
c_p = 920 # J/kgK; specific heat
T_0 = T_air # oC
I_0 = 500 # W/m^2; insolation value
alpha = 0.95 # fraction of radiation absorbed
C = rho * c_p * A * (L / n) # J/oC; capacitance
R_k = (1 / (k * A)) * (L / n) # oC/W; resistance
Q_dot_solar = alpha * I_0 * A # W; absorbtion rate
epsilon = 0.9 # emissivity
sigma = 5.669e-8 # W/m^2K^4; Stefan-Boltzmann constant
h_rad = epsilon * sigma * ((T_0 + 273)^2 + (T_sky + 273)^2) * ((T_0 + 273) + (T_sky + 273)) # W/m^2oC; radiative coefficient
h_conv = 7 # W/m^2K; convection coefficient
R_rad = 1 / (h_rad * A) # oC/W; radiation resistance
R_conv = 1 / (h_conv * A) # oC/W; convection resistance

dt_max_surface = R_k * C / 3 # s; maximum allowable step size for the surface node
dt_max_interior = R_k * C / 2 # s; maximum allowable step size for inteior nodes
dt_max_bottom = R_k * C # s; maximum allowable step size for bottom node

'''SIMULATION'''
# Setup
t = np.arange(0, 60 * 60 * 2, dt) # s
p = 0
t_p = t[0] # s; current time
T = np.zeros((n + 1, len(t))) # oC
T[:, 0] = T_0

# Loop
while (t_p < t[-1]):
    # Surface
    # Sunny (1st hour)
    if (t_p <= 60 * 60):
        T[0, p] = ((alpha * I_0 * A) + (T_sky / R_rad) + (T_air / R_conv) + (T[1, p] / (R_k / 2))) / ((1 / R_rad) + (1 / R_conv) + (1 / (R_k / 2))) 
    # Cloudy (2nd hour)
    elif (t_p > 60 * 60):
        T[0, p] = ((T_cloud / R_rad) + (T_air / R_conv) + (T[1, p] / (R_k / 2))) / ((1 / R_rad) + (1 / R_conv) + (1 / (R_k / 2)))
    
    # Top Node
    T[1, p + 1] = T[1, p] + (dt / C) * (((T[2, p] - T[1, p]) / R_k) + ((T[0, p] - T[1, p]) / (R_k / 2))) 
    
    # Interior Node(s)
    i = 2
    while (i <= n - 1):
        T[i, p + 1] = T[i, p] + (dt / C) * (((T[i - 1, p] - T[i, p]) / R_k) + ((T[i + 1, p] - T[i, p]) / R_k))
        
        i = i + 1
    
    # Bottom Node
    T[n, p + 1] = T[n, p] + (dt / C) * ((T[n - 1, p] - T[n, p]) / R_k)
    
    p = p + 1
    t_p = t[p]

'''PLOTS'''
# Temperature vs Time
plt.figure()
for i in range(n + 1):
    if (i == 0):
        plt.plot(t, T[i, :], label = 'Surface')
    else:
        plt.plot(t, T[i, :], label = 'Node: ' + str(i))
plt.title('Temperature vs Time')
plt.xlabel('Time [Seconds]')
plt.xlim(0, 60 * 60 * 2 - dt)
plt.ylabel('Temperature [Degrees-Celcius]')
plt.ylim(30, 60)
plt.legend()
plt.show()

# Temperature vs Position
x_0 = np.array([0])
x = np.arange(L / n / 2, L - (L / n / 2), L / n) # m
x = np.concatenate((x_0, x), axis = 0)
i = int(60 * 60 * 2 / 10)
T = np.column_stack((T[:, 0 * i], T[:, 1 * i], T[:, 2 * i], T[:, 3 * i], T[:, 4 * i], T[:, 5 * i], T[:, 6 * i], T[:, 7 * i], T[:, 8 * i], T[:, 9 * i], T[:, 10 * i]))
plt.figure()
for j in range(len(x)):
    plt.plot(x, T[:, j], label = 't = ' + str(j * 10))
plt.title('Temperature vs Position')
plt.xlabel('Position [Meters]')
plt.ylabel('Temperature [Degrees-Celcius]')
plt.legend()














