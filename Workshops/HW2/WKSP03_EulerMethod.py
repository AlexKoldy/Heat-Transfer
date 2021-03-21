# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 12:06:27 2021

@author: george.sidebotham2
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# SETUP EXCEL -> PYTHON
data = pd.read_excel (r'CoffeeData.xlsx', sheet_name=('DATA'))
#imports data from excel file
data.to_csv('CoffeeData.txt')
#saves data in a csv file
data = np.loadtxt('CoffeeData.txt', dtype = np.float, delimiter = ',', skiprows=1)
#imports data to a list
#  Convert data from list to individual arrays (for data analysis)    
time = data[:,1]    #note, 1st column (index 0) is a  list of numbers
T_top = data[:,2]
T_side = data[:,3]

# CONSTANTS
tau = 60 # Goal Seek Tau: 49.5770152834528 # min

# INITIAL PARAMETERS
T_0 = 75 #Degrees Celcius
T_inf = 21 #Degrees Celcius
C = 2010 #J/K
A = 0.0496 #m^2
U = C / (tau * 60) / A # W/(m^2-K)
dt = 0.5 # min
t0 = 0 # min
t = np.zeros(501)
T = np.zeros(501)
T[0] = T_0

# SIMULATION
p = 0
while (p <= t.shape[0] - 1):
    if (p != t.shape[0] - 1):
        t[p + 1] = t[p] + dt
        T[p + 1] = T[p] + dt * ((T_inf - T[p]) / (tau)) #Euler Method
        
    p = p + 1

# PLOTS
plt.figure(1)
plt.plot(time,T_top,'o', label = 'Top')
plt.plot(time,T_side, '+', label = 'Side')
plt.plot(t, T, label = 'Euler Method')
plt.xlabel('time (min)'), plt.ylabel('T_top (oC)')
plt.legend()
plt.title('Mug Temperature Simulation Curve and Data (Tau = 60)')