'''Heat Transfer Workshop 11
Alexander Koldy, William Lin, Solbi Oh '''

'''Constants'''
# Air
k_air = 0.0257 # W/mK; themal conductivity
rho_air = 1.205 # kg/m^3; density
c_p_air = 1009 # J/kgK; specific heat
mu_air = 1.82 * 10**(-5) # kg/ms; dynamic viscocity
beta_air = 3.41 * 10**(-3) # 1/K; coefficient of thermal expansion

# Water
k_H2O = 0.6 # W/mK; themal conductivity
rho_H2O = 1000 # kg/m^3; density
c_p_H2O = 4200 # J/kgK; specific heat
mu_H2O = 1 * 10**(-3) # kg/ms; dynamic viscocity
beta_H2O = 2.07 * 10**(-4) # 1/K; coefficient of thermal expansion

'''Input Parameters'''
D = float(input("Please input an outer diameter (D) [m]: ")) # m; outer diameter
V = float(input("Please input a wind speed (V) [m/s]: ")) # m/s; wind speed
T_inf = float(input("Pleae input an ambient temperature (T_inf) [oC]: ")) # oC; ambient temperature
fluid = input("Please input the type of fluid ('air' or 'water'): ") # type of fluid

'''Setup'''
if (fluid == 'air'):
    k = k_air
    rho = rho_air
    c_p = c_p_air
    mu = mu_air
    beta = beta_air
elif (fluid == 'water'):
    k = k_H2O
    rho = rho_H2O
    c_p = c_p_H2O
    mu = mu_H2O
    beta = beta_H2O
    
T_surface = 32 # oC; surface temperature of skin (only for natural convection)
T_film = (T_surface + T_inf) / 2 # oC; film temperature

g = 9.8 # m/s; acceleration due to gravity
v = mu / rho # m^2/s
Re = rho * V * D / mu # Reynolds number
alpha = k / (rho*c_p) # m^2/s; thermal diffusivity
Ra = (g * beta * abs(T_surface - T_inf) * D**3) / (v * alpha) # Rayleigh Number 


'''Hilpert/Knudesn/Katz Correlation'''
def Hilpert_Knudesn_Katz():
    # Find C and n by using Reynolds number range
    if (0.4 <= Re and Re <= 4):
        C = 0.989
        n = 0.330
    elif (4 <= Re and Re <= 35):
        C = 0.911
        n = 0.385
    elif (35 <= Re and Re <= 4083):
        C = 0.683
        n = 0.466
    elif (4083 <= Re and Re <= 40045):
        C = 0.193
        n = 0.618
    elif (40045 <= Re and Re <= 400000):
        C = 0.0266
        n = 0.805
    
    return C * (rho*V*D / mu)**n * (v / alpha)**(1/3) * k / D # W/m^2K; convection coefficient

'''Notsurewho Correlation'''
def Notsurewho():
    # Find C and N by using Rayleigh number range
    if (10**4 <= Ra and Ra <= 2.12 * 10**7):
        C = 0.53
        N = 0.25
    elif (2.12 * 10**7 <= Ra and Ra <= 10**12):
        C = 0.13
        N = 0.3333
        
    return C * Ra**N * k / D # W/m^2K; convection coefficient

'''Final Output'''
if (V == 0):
    h = Notsurewho()
else:
    h = Hilpert_Knudesn_Katz()
    
print("Convection Coefficient (h): ")
print(h)

