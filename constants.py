# importing the modules
import numpy as np
import matplotlib.pyplot as plt

### Physcial Constants in cm
kbT = 0.025852 
q = 1.602117*np.power(10.0, -19)
e_O = 8.854187*np.power(10.0, -14) #cm
e_s = 11.7*e_O
e_ox = 3.9*e_O
root_2 = np.power(2.0, 0.5)

### Basic Semiconductor Constants
n_i = 1.07*np.power(10.0, 10)
n_i_2 = np.power(n_i, 2)
N_v = 3.1*np.power(10.0, 19)
N_c = 2.86*np.power(10.0, 19)
X_si = 4.05
Eg = 1.1242
B = 1/kbT