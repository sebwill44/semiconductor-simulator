# importing the modules
import numpy as np
import matplotlib.pyplot as plt
import printing as p
import recombo
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

u_p_P = lambda Na : 49.7 + (418.3/(1+np.power(Na/(1.6*np.power(10.0, 17)), .7)))
u_n_P = lambda Na : 232 + (1180/(1+np.power(Na/(8.0*np.power(10.0, 16)), .9)))
u_n_N = lambda Nd : 92 + (1268/(1+np.power(Nd/(1.3*np.power(10.0, 17)), .91)))
u_p_N = lambda Nd : 130 + (370/(1+np.power(Nd/(8.0*np.power(10.0, 17)), 1.25)))

D_n_P = lambda Na : kbT*u_n_P(Na)
D_p_N = lambda Nd : kbT*u_p_N(Nd)

L_p_N = lambda Nd : np.sqrt(recombo.T_rec_N(Nd)*D_p_N(Nd))
L_n_P = lambda Na : np.sqrt(recombo.T_rec_P(Na)*D_n_P(Na))

""" x = np.linspace(np.power(10.0, 14), np.power(10.0, 20))
plt.plot(x, u_n_P(x), 'g')
plt.plot(x, u_p_P(x),'g--')
plt.plot(x, u_n_N(x), 'r')
plt.plot(x, u_p_N(x),'r--')
plt.xscale('log')
plt.title("Minority-carrier Recombination Lifetime")
plt.ylabel("Carrier Recombination Lifetime (s)")
plt.xlabel("Position $x$ ($cm$) ")
plt.grid()
plt.legend([r'$\mu_{n,P}(N_{a}^{-})$', r'$\mu_{p,P}(N_{a}^{-})$', r'$\mu_{n,N}(N_{d}^{+})$', r'$\mu_{p,N}(N_{d}^{+})$'])

plt.show() """
