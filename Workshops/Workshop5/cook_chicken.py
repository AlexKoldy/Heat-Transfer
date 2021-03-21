from chicken import Chicken
from methods import Methods

T_1_0 = 5 # oC
T_2_0 = 5 # oC
f = 0.001

boiling = Methods('Boiling')
baking = Methods('Baking')
bbqing = Methods('BBQing')

chicken_1 = Chicken()
chicken_1.prepare(T_1_0, T_2_0, boiling.T_inf, boiling.h, f)
chicken_1.cook('Boiling')

chicken_2 = Chicken()
chicken_2.prepare(T_1_0, T_2_0, baking.T_inf, baking.h, f)
chicken_2.cook('Baking')

chicken_3 = Chicken()
chicken_3.prepare(T_1_0, T_2_0, bbqing.T_inf, bbqing.h, f)
chicken_3.cook('BBQing')

