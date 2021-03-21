class Methods:
    def __init__(self, method):        
        if (method == 'Boiling'): 
            self.h = 1000 # W/m^2oC
            self.T_inf = 100 # oC
            
        elif (method == 'Baking'): 
            self.h = 15 # W/m^2oC
            self.T_inf = 177 # oC
            
        elif (method == 'BBQing'): 
            self.h = 25  # W/m^2oC
            self.T_inf = 500 # oC
        
