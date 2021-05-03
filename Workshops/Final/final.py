'''Heat Transfer Final Project: 
Simulating Heat Transfer Inside Mycelium Container
Alexander Koldy, William Lin, Solbi Oh'''
import numpy as np
import matplotlib.pyplot as plt

'''===========CONSTANTS=========='''
sigma = 5.669e-8 # W/m^2K^4; Stefan-Boltzmann Constant
epsilon= 0.92 # emissivity; (Reference: https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html)

'''----------Dimensions & Areas----------
in: inner area
out: outer area'''
# Glass Cup
h_glass = 0.0381 # m; height of glass including bottom thickness
d_glass = 0.0508 # m; outer diameter of glass
t_glass_side = 0.004 # m
t_glass_bot = 0.015 # m
A_glass_side_in = np.pi * (d_glass - 2*t_glass_side) * (h_glass - t_glass_bot) # m^2
A_glass_side_out = np.pi * d_glass * h_glass # m^2
A_glass_bot_in = (np.pi/4) * (d_glass - 2*t_glass_side)**2 # m^2
A_glass_bot_out = (np.pi/4) * d_glass**2 # m^2

# Mycelium Box
h_myc = 0.12065 # m
w_myc = 0.1397 # m
l_myc = 0.12065 # m
t_myc = 0.01905 # m
A_myc_top_in = (w_myc - 2*t_myc) * (l_myc - 2*t_myc) # m^2
A_myc_top_out = w_myc * l_myc # m^2
A_myc_side_in = 2 * (w_myc - 2*t_myc) * (h_myc - 2*t_myc) + 2 * (l_myc - 2*t_myc) * (h_myc - 2*t_myc) # m^2
A_myc_side_out = 2* (w_myc * h_myc) + 2 * (l_myc * h_myc) # m^2
A_myc_bot_in = (w_myc - 2*t_myc) * (l_myc - 2*t_myc) # m^2
A_myc_bot_out = w_myc * l_myc # m^2

'''----------Heat & Density Properties----------'''
# Water
h_water = 50 # W/m^2K
h_water_air = 50 # W/m^2K
c_p_water = 4184 # J/kgK; (Reference: https://en.wikipedia.org/wiki/Specific_heat_capacity#:~:text=The%20SI%20unit%20of%20specific,%E2%88%921%E2%8B%85K%E2%88%921)
rho_water = 997 # kg/m^3; (Reference: https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html)

# Glass Cup
h_water_glass = 50 # W/m^2K
h_air_glass = 5 # W/m^2K
k_glass = 0.8 # W/mK
c_p_glass = 840 # J/kgK; (Reference: https://www.engineeringtoolbox.com/specific-heat-solids-d_154.html)
rho_glass = 2500 # kg/m^3; (Reference: https://uk.saint-gobain-building-glass.com/en-gb/architects/physical-properties#:~:text=The%20density%20of%20glass%20is,or%202500%20kg%20per%20m3.&text=The%20compressive%20strength%20of%20glass,N%2Fmm2%20%3D%201000%20MPa)

# Air
c_p_air = 1005 # J/kgK; (Reference: https://www.engineeringtoolbox.com/air-specific-heat-capacity-d_705.html)
rho_air = 1.225 # kg/m^3; (Reference: https://www.macinstruments.com/blog/what-is-the-density-of-air-at-stp/#:~:text=According%20to%20the%20International%20Standard,%3A%200.0765%20lb%2Fft%5E3)

# Mycelium Box
h_air_myc = 5 # W/m^2K; (Reference: http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thrcn.html)
k_myc = 0.05 # W/mK; (Reference: https://www.biorxiv.org/content/10.1101/569749v1.full.pdf)
h_myc_atm = 5 # W/m^2K
c_p_myc = 10200 # J/kgK; (Reference: "Physico-Mechanical and Thermodynamic Properties of Mycelium-Based Biocomposites: A Review")
rho_myc = 170 # kg/m^3; (Reference: https://www.sciencedirect.com/science/article/pii/S0264127518308347#:~:text=3.2.&text=Density%20of%20the%20mycelium%20composite,pressed%20mycelium%20composites%20%5B7%5D)

'''----------Resistances----------
Relative to middle node in model
R_h: R_h_<MEDIUM_1>_<MEDIUM_2>_<bot/side/top>
R_k: R_k_<MEDIUM>_<bot/side/top>
R_r: R_r_<MEDIUM_1>_<MEDIUM_2>_<bot/side/top>'''
# Bottom
R_h_water_glass_bot = 1 / (h_water_glass * A_glass_bot_in) # oC/W
R_k_glass_bot = t_glass_bot / (k_glass * A_glass_bot_out) # oC/W
R_k_myc_bot = t_myc / (k_myc * A_myc_bot_out) # oC/W
R_h_myc_atm_bot = 1 / (h_myc_atm * A_myc_bot_out) # oC/W

# Side
R_h_water_glass_side = 1 / (h_water_glass * A_glass_side_in) # oC/W
R_k_glass_side = t_glass_side / (k_glass * A_glass_side_in) # oC/W
R_h_glass_air_side = 1 / (h_air_glass * A_glass_side_out) # oC/W
R_h_air_myc_side = 1 / (h_air_myc * A_myc_side_in) # oC/W
R_k_myc_side = t_myc / (k_myc * A_myc_side_out) # oC/W
R_h_myc_atm_side = 1 / (h_myc_atm * A_myc_side_out) # oC/W

# Top
R_h_water_top = 1 / (h_water * A_glass_bot_in) # oC/W
R_h_water_air_top = 1 / (h_water_air * A_glass_bot_in) # oC/W
R_h_air_myc_top = 1 / (h_air_myc * A_myc_top_in) # oC/W
R_k_myc_top = t_myc / (k_myc * A_myc_top_out) # oC/W
R_h_myc_atm_top = 1 / (h_myc_atm * A_myc_top_out) # oC/W

'''==========SIMULATIONS=========='''
#def 1_node_model():
    #continue

def six_node_model():
    '''----------Setup----------'''
    steps = 5000001 # total amount of timesteps
    t = np.zeros(steps) # s
    dt = 0.01 # s
    t_sim = steps * dt # s; total simulation time
    p = 0
    
    '''----------Capacitance----------'''
    # Water Node
    C_water = c_p_water * rho_water * A_glass_bot_in * (h_glass - t_glass_bot) # J/K
    
    # Glass Cup Nodes
    C_glass_bot = c_p_glass * rho_glass * A_glass_bot_out * t_glass_bot # J/K
    C_glass_side = c_p_glass * rho_glass * (A_glass_bot_out - A_glass_bot_in) * (h_glass - t_glass_bot) # J/K
    
    # Air
    C_air_side = c_p_air * rho_air * (A_myc_bot_in - A_glass_bot_out) * (h_glass) # J/K
    C_air_top = c_p_air * rho_air * A_myc_bot_in * (h_myc - 2*t_myc - h_glass) # J/K
    
    # Mycelium Nodes
    C_myc_bot = c_p_myc * rho_myc * A_myc_bot_out * t_myc # J/K
    C_myc_side = c_p_myc * rho_myc * (A_myc_bot_out - A_myc_bot_in) * (h_myc - 2*t_myc) # J/K
    C_myc_top = c_p_myc * rho_myc * A_myc_top_out * t_myc # J/K
    
    '''----------Temperatures----------'''
    # Water Node
    T_water = np.zeros(steps) # oC
    T_water[0] = 32.2 # oC
    
    # Glass Cup Nodes
    T_glass_bot = np.zeros(steps) # oC
    T_glass_side = np.zeros(steps) # oC
    T_glass_bot[0] = 20 # oC
    T_glass_side[0] = 20 # oC

    # Myeclium Nodes
    T_myc_bot = np.zeros(steps) # oC
    T_myc_side = np.zeros(steps) # oC
    T_myc_top = np.zeros(steps) # oC
    T_myc_bot[0] = 20 # oC
    T_myc_side[0] = 20 # oC
    T_myc_top[0] = 20 # oC
    
    # Enviornment
    T_inf = 20 # oC; ambient temperature
    
    '''----------Simulation----------'''
    while p < steps - 1:
        # Update radiation
        h_r_glass_myc = epsilon * sigma * (T_glass_side[p]**2 + T_myc_side[p]**2) * (T_glass_side[p] + T_myc_side[p]) # W/m^2K
        R_r_glass_myc_side = 1 / (h_r_glass_myc * A_glass_side_out) # oC/W

        # Main node
        T_water[p + 1] = T_water[p] + (dt/C_water) * (((T_glass_bot[p] - T_water[p]) / (R_h_water_glass_bot + R_k_glass_bot/2)) + ((T_glass_side[p] - T_water[p]) / (R_h_water_glass_side + R_k_glass_side/2)) + ((T_myc_top[p] - T_water[p]) / (R_h_water_top + R_h_water_air_top + R_h_air_myc_top + R_k_myc_top/2)))
        
        # Bottom
        T_glass_bot[p + 1] = T_glass_bot[p] + (dt/C_glass_bot) * (((T_water[p] - T_glass_bot[p]) / (R_h_water_glass_bot + R_k_glass_bot/2)) + ((T_myc_bot[p] - T_glass_bot[p]) / (R_k_glass_bot/2 + R_k_myc_bot/2)))
        T_myc_bot[p + 1] = T_myc_bot[p] + (dt/C_myc_bot) * (((T_glass_bot[p] - T_myc_bot[p]) / (R_k_glass_bot/2 + R_k_myc_bot/2)) + ((T_inf - T_myc_bot[p]) / (R_k_myc_bot/2 + R_h_myc_atm_bot)))
        
        # Side
        T_glass_side[p + 1] = T_glass_side[p] + (dt/(C_glass_side + C_air_side + C_air_top)) * (((T_water[p] - T_glass_side[p]) / (R_h_water_glass_side + R_k_glass_side/2)) + ((T_myc_side[p] - T_glass_side[p]) / (R_k_glass_side/2 + (1 / ((1 / (R_h_glass_air_side + R_h_air_myc_side)) + (1 / R_r_glass_myc_side))) + R_k_myc_side/2)))
        T_myc_side[p + 1] = T_myc_side[p] + (dt/C_myc_side) * (((T_glass_side[p] - T_myc_side[p]) / (R_k_glass_side/2 + (1 / ((1 / (R_h_glass_air_side + R_h_air_myc_side)) + (1 / R_r_glass_myc_side))) + R_k_myc_side/2)) + ((T_inf - T_myc_side[p]) / (R_k_myc_side/2 + R_h_myc_atm_side)))
    
        # Top
        T_myc_top[p + 1] = T_myc_top[p] + (dt/C_myc_top) * (((T_water[p] - T_myc_top[p]) / (R_h_water_top + R_h_water_air_top + R_h_air_myc_top + R_k_myc_top/2)) + ((T_inf - T_myc_top[p]) / (R_k_myc_top/2 + R_h_myc_atm_top)))
        
        # Increment
        t[p + 1] = t[p] + dt
        p += 1
        
        # do we need to consider entire capacitance or only the indivial note capacitance?
        # should we lump both top and side air to the side of the mug?
        # what values to put for convective coefs? (there is a large range)
        # fluids in series? 
    plt.figure()
    plt.plot(t, T_water)

six_node_model()
    


 