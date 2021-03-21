class Methods:
    def __init__(self, method):
        self.T_1_0 = 5 # oC
        self.T_2_0 = 5 # oC
        self.f = 0.001
        
        if (method == 'Boiling'): 
            self.h = 1000 #15, 25]  # W/m^2oC
            self.T_inf = 100 # 177, 500] # oC
            
        if (method == 'Baking'): 
            self.h = 15 #25]  # W/m^2oC
            self.T_inf = 177 #, 50 # oC
        if (method == 'BBQing'): 
            self.h = 25]  # W/m^2oC
            self.T_inf = 500 # oC
