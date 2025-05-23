''' 
    exercise_3_1.py
    
    Author: Aldo Nunez Tovar
'''

import semiconductor as sc
from astropy import units as u 

T_1 = 50 * u.Kelvin
T_2 = 350 * u.Kelvin

n_i_1 = sc.n_i(T_1)
n_i_2 = sc.n_i(T_2)

print(f"Exercise 3.1:")
print(f"n_i({T_1:.0f}) = {n_i_1:0.2E}")
print(f"n_i({T_2:.0f}) = {n_i_2:0.2E}")
