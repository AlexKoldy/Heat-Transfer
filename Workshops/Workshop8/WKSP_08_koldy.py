import matplotlib.pyplot as plt
import numpy as np

'''STEP 1: FIN DESIGN'''
# Constants and Inputs
efficiency_fin = 0.75 # heat fin efficency
k = 200 # W/mK; conductive coefficent of aluminum 
h_air = 25 # W/m^2K; convective coefficent of air
w_fin = 0.1 # m; width of the fin
d_tube = 0.03 # m; diameter of pipe

# Setup
t_fin = np.linspace(0.0001, 0.001, 1000) # m; thickness of fin
L_fin = (w_fin / 2) - (d_tube / 2) # m; length of fin 
P_fin = 2 * (t_fin + w_fin) # m; perimeter of fin
A_fin = t_fin * w_fin # m^2; cross-secitonal area of fin
m = (h_air * P_fin / (k * A_fin))**0.5
P_avg_fin = (d_tube * np.pi + 4 * w_fin) * 0.5 # m; average circumference/perimeter
A_avg_fin = P_avg_fin * t_fin # m^2; average croessectional area
L_avg_fin = (0.5 * (0.5 + (1 / (2**0.5)))) * w_fin - (d_tube / 2) # m; average distance between edge of fin and pipe
efficiency = 1 / (((m * L_fin)**2 / 2) + 1)

# Plots
plt.figure()
plt.title('Efficency vs Fin Thickness')
plt.xlabel('Thickness [m]')
plt.ylabel('Efficiency')
plt.plot(t_fin, efficiency)
plt.axhline(y = efficiency_fin, color = 'black', label = 'Heat Fin')
plt.legend()

'''STEP 2: OVERALL HEAT TRANSFER COEFFICIENT'''
# Constants and Inputs
t_tube = 0.003 # m; thickness of the pipe
t_fin = 0.0005 # m; thickness
k = 200 # W/mK; conductive coefficent of aluminum
dx = 0.01 # m; spacing between fins
h_oil = 200 # W/m^2K; convective coefficent of oil
h_air = 25 # W/m^2K; convective coefficent of air
w_fin = 0.1 # m; width of the fin
d_tube = 0.03 # m; diameter of pipe

# Setup
A_inner_tube = 2 * np.pi * (d_tube / 2 - t_tube) * dx # m^2; inner surface area of pipe
A_outer_tube = np.pi * d_tube * dx # m^2; outer surface area of pipe
A_avg_tube = (A_outer_tube + A_inner_tube) / 2 # m^2; average surface area of pipe
A_tube_fin = np.pi * d_tube * t_fin # m^2; area between tube and fin
L_fin = (w_fin / 2) - (d_tube / 2) # m; length of fin
A_surface_fin = 2 * (w_fin**2) + 4 * w_fin * t_fin # m^2; surface area of fin
P_avg_fin = (d_tube * np.pi + 4 * w_fin) * 0.5 # m; average circumference/perimeter
A_avg_fin = P_avg_fin * t_fin # m^2; average croessectional area
L_avg_fin = (0.5 * (0.5 + (1 / (2**0.5)))) * w_fin - (d_tube / 2) # m; average distance between edge of fin and pipe
R_h_oil = 1 / (h_oil * A_inner_tube) # oC/W; convective resistance through oil
R_k_tube = t_tube / (k * A_avg_tube) # oC/W; conductive resistance through pipe
R_h_air_tube = 1 / (h_air * (A_outer_tube - A_tube_fin)) # oC/W; convective resistance between air and pipe
R_k_fin = (L_avg_fin / 2) / (k * A_avg_fin) # oC/W; conductive resistance through fin
R_h_air_fin = 1 / (h_air * A_surface_fin) # oC/W; convective resistance between air and fin
R_total = R_h_oil + R_k_tube + (1 / ((1 / (R_k_fin + R_h_air_fin)) + (1 / R_h_air_tube))) # oC/W; equivalent resistance

U = 1 / (R_total * A_inner_tube)

# Output
print('R_h_oil: ' + str(R_h_oil) + ' oC/W')
print('R_k_tube: ' + str(R_k_tube) + ' oC/W')
print('R_h_air_tube: ' + str(R_h_air_tube) + ' oC/W')
print('R_k_fin: ' + str(R_k_fin) + ' oC/W')
print('R_h_air_fin: ' + str(R_h_air_fin) + ' oC/W')
print('R_total: ' + str(R_total) + ' oC/W')
print('U: ' + str(U) + ' W/m^2K')


