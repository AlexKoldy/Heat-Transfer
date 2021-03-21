import numpy as np
import matplotlib.pyplot as plt

class Chicken:
    def __init__(self):
        self.rho = 1100  # kg/m^3
        self.c_p = 3200 # J/kgoC
        self.L = 0.019 # m
        self.A = 0.002 # m^2
        self.k = 0.46 # W/mK
        
        self.T_cooked = 71.40 # oC
        
        self.is_cooked_1 = False
        self.is_cooked_2 = False
        self.is_cooked_center = False
        self.is_cooked = self.is_cooked_center
                
    def prepare(self, T_1_0, T_2_0, T_inf, h, f):
        self.T_1_0 = T_1_0 # oC
        self.T_2_0 = T_2_0 # oC
        self.T_inf = T_inf # oC
        
        self.t = np.zeros(50001) # s
        self.T_1 = np.zeros(50001) # oC
        self.T_2 = np.zeros(50001) # oC
        self.T_center = np.zeros(50001) # oC
        self.T_surface = np.zeros(50001) # oC
        self.T_1[0] = self.T_1_0;
        self.T_2[0] = self.T_2_0;
        
        self.C_1 = self.rho * self.c_p * (self.L / 2) * self.A # J/oC
        self.C_2 = self.rho * self.c_p * (self.L / 2) * self.A # J/oC
        self.R_12 = (self.L / 2) / (self.k * self.A) # oC/W
        self.R_2s = (self.L / 4) / (self.k * self.A) # oC/W
        self.R_h = 1 / (h * self.A) #oC/W
        
        self.tau_1 = self.C_1 * self.R_12 # s
        self.tau_2 = self.C_2 / ((1 / self.R_12) + (1 / (self.R_2s + self.R_h))) # s
        
        self.dt =  f * min(self.tau_1, self.tau_2) # s

    def cook(self, method):
        print("Cooking chicken using: " + method)
        
        if (method == 'Boiling'): 
            dT = 0.002
            
        elif (method == 'Baking'): 
            dT = 0.005
            
        elif (method == 'BBQing'): 
            dT = 0.05
            
        p = 0
        is_done = False
                    
        while (p <= self.t.shape[0] - 1):

            
            if (self.is_cooked == True and is_done == False):
                self.t_end = self.t[p]
                is_done = True
            
            if (p != self.t.shape[0] - 1):
                self.t[p + 1] = self.t[p] + self.dt
        
                # Euler Method
                self.T_1[p + 1] = self.T_1[p] + ((self.dt / self.C_1) * ((self.T_2[p] - self.T_1[p]) / self.R_12))
                self.T_2[p + 1] = self.T_2[p] + ((self.dt / self.C_1) * (((self.T_1[p] - self.T_2[p]) / self.R_12) + ((self.T_inf - self.T_2[p]) / (self.R_2s + self.R_h)))) 
      
            # Surface Temperature Calculations
            self.T_surface[p] = self.T_inf + ((self.R_h / (self.R_2s + self.R_h)) * (self.T_2[p] - self.T_inf))
          
            # Center Temperature Calculations
            self.T_center[p] = ((9 * self.T_1[p]) - self.T_2[p]) / 8
        
            # Intersction Estimate
            self.check_cooked(self.T_1[p], '1', self.t[p], self.is_cooked_1, dT, p)
            self.check_cooked(self.T_2[p], '2', self.t[p], self.is_cooked_2, dT, p)
            self.check_cooked(self.T_center[p], 'center', self.t[p], self.is_cooked_center, dT, p)
            
            p = p + 1
            
        self.plot(method)

    def check_cooked(self, T, location, t, is_cooked, dT, p):
        if (T < self.T_cooked + dT and T > self.T_cooked - dT):
            print("Location = ", location, " cooked at: t = ", t / 60, " minutes")
            
            if (location == 'center'):
                print("Surface temperature (degrees-Celcius) when center is cooked: ", self.T_surface[p])
                self.is_cooked = True
                

            
    def plot(self, method):
        plt.figure()
        plt.title('Temperature of Chicken using ' + method + ' [Temperature vs Time]')
        plt.xlabel('Time [minutes]')
        plt.ylabel('Temperature [degrees-Celcius]')
        plt.xlim(0, 200)
        plt.ylim(0, 500)

        # Convert time to minutes when plotting
        plt.plot(self.t / 60, self.T_1)
        plt.plot(self.t / 60, self.T_2)
        plt.plot(self.t / 60, self.T_surface)
        plt.plot(self.t / 60, self.T_center)
        plt.axhline(self.T_cooked, color='black')
        plt.legend(['T_1', 'T_2', 'T_surface', 'T_center', 'T_Cook'])
        plt.show()
        
        # Setup Temperature vs. Position by using increments
        i = 100
        t = [self.t[i], self.t[2 * i], self.t[3 * i], self.t[4 * i], self.t[5 * i], self.t[6 * i], self.t[7 * i], self.t[8 * i], self.t[9 * i], self.t[10 * i]]
        T_1 = [self.T_1[i], self.T_1[2 * i], self.T_1[3 * i], self.T_1[4 * i], self.T_1[5 * i], self.T_1[6 * i], self.T_1[7 * i], self.T_1[8 * i], self.T_1[9 * i], self.T_1[10 * i]]
        T_2 = [self.T_2[i], self.T_2[2 * i], self.T_2[3 * i], self.T_2[4 * i], self.T_2[5 * i], self.T_2[6 * i], self.T_2[7 * i], self.T_2[8 * i], self.T_2[9 * i], self.T_2[10 * i]]
        T_center = [self.T_center[i], self.T_center[2 * i], self.T_center[3 * i], self.T_center[4 * i], self.T_center[5 * i], self.T_center[6 * i], self.T_center[7 * i], self.T_center[8 * i], self.T_center[9 * i], self.T_center[10 * i]]
        T_surface = [self.T_surface[i], self.T_surface[2 * i], self.T_surface[3 * i], self.T_surface[4 * i], self.T_surface[5 * i], self.T_surface[6 * i], self.T_surface[7 * i], self.T_surface[8 * i], self.T_surface[9 * i], self.T_surface[10 * i]]
        pos = [0, self.L / 4, 3 * self.L / 4, self.L]
        
        plt.figure()
        plt.title('Temperature of Chicken using ' + method + ' [Temperature vs Position]')
        plt.xlabel('Position [m]')
        plt.ylabel('Temperature [degrees-Celcius]')
        
        for j in range(10):
            plt.plot(pos, [T_center[j], T_1[j], T_2[j], T_surface[j]], label=(str(round(t[j] / 60, 2)) + ' mins'))
            
        print(T_center[j])
        plt.legend()  
        plt.show()
        
        


