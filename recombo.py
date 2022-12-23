# importing the modules
import numpy as np
import matplotlib.pyplot as plt
import printing as p

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

T_rec_P = lambda Na : 1/((3.45*np.power(10.0, -12)*Na)+(9.5*np.power(10.0, -32)*Na*Na))
T_rec_N = lambda Nd : 1/((7.8*np.power(10.0, -13)*Nd)+(1.8*np.power(10.0, -31)*Nd*Nd))
T_gen_p = lambda Na : 75*T_rec_P(Na)
T_gen_n = lambda Nd : 75*T_rec_N(Nd)

""" x = np.linspace(np.power(10.0, 14), np.power(10.0, 20))
plt.plot(x, T_rec_N(x), 'g')
plt.plot(x, T_rec_P(x),'r')
plt.yscale('log')
plt.xscale('log')
plt.title("Minority-carrier Recombination Lifetime")
plt.ylabel("Carrier Recombination Lifetime (s)")
plt.xlabel("Position $x$ ($cm$) ")
plt.grid()
plt.legend([r'$\tau_{rec,N}(N_{d}^{+})$', r'$\tau_{rec,P}(N_{a}^{+})$'])

plt.show() """
