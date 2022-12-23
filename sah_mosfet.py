# importing the modules
import numpy as np
import matplotlib.pyplot as plt

### Physcial Constants
kbT = 0.025852 #eV
q = 1.602*np.power(10.0, -19)
epsilon_s = 11.7*(8.84*np.power(10.0,-14))
epsilon_sio2 = 3.9*(8.84*np.power(10.0,-14))
root_2 = np.power(2.0, 0.5)

### Basic Semiconductor Constants
n_i = 1.05*np.power(10.0, 10)
n_i_2 = np.power(n_i, 2)
N_v = 3.1016*np.power(10.0, 19)
X_si = 4.05
Eg = 1.1242
B = 1/kbT

t_ox = 4*np.power(10.0, -7)
C_ox = 8.632*np.power(10.0, -7)
L_ch = 0.2*np.power(10.0, -4)
W_ch = 15*np.power(10.0, -4)

p_t_ox = 4*np.power(10.0, -7)
p_C_ox = 8.632*np.power(10.0, -7)
p_L_ch = 0.2*np.power(10.0, -4)
p_W_ch = 37.5*np.power(10.0, -4.)
Vtn = 1.5

def sah_I_D_lin_nmos(Vds, Vgs):
    return 500*C_ox*(W_ch/L_ch)*((Vgs-Vtn)*Vds-0.5*Vds*Vds)

def sah_I_D_sat_nmos(Vds, Vgs):
    return 500*C_ox*(W_ch/(2*L_ch))*(np.power(Vgs-Vtn, 2))

def sah_I_D_lin_pmos(Vds, Vgs):
    return 200*p_C_ox*(p_W_ch/p_L_ch)*((Vgs-Vtn)*Vds-0.5*Vds*Vds)

def sah_I_D_sat_pmos(Vds, Vgs):
    return 200*p_C_ox*(p_W_ch/(2*p_L_ch))*(np.power(Vgs-Vtn, 2))

### If velocity saturation is considered

ne_sat = 2*(4*np.power(10.0, 6))/500
pe_sat = 2*(4*np.power(10.0, 6))/200


def V_dsat(Vgs):
    return Vgs-Vtn

e_sat = 2*(4*np.power(10.0, 6))/500

def vsat_sah_I_D_sat_nmos(Vds, Vgs):
    return 4*np.power(10.0, 6)*C_ox*W_ch*(np.power(Vgs-Vtn, 2)/((Vgs-Vtn)+ne_sat*L_ch))

def vsat_sah_I_D_sat_pmos(Vds, Vgs):
        return 4*np.power(10.0, 6)*C_ox*p_W_ch*(np.power(Vgs-Vtn, 2)/((Vgs-Vtn)+pe_sat*L_ch))

def nsat_Vdsat(Vgs):
    return ((Vgs-Vtn)*ne_sat*L_ch)/((Vgs-Vtn)+L_ch*ne_sat)

def psat_Vdsat(Vgs):
    return ((Vgs-Vtn)*pe_sat*L_ch)/((Vgs-Vtn)+L_ch*pe_sat)

print("V_dsat: " + str(V_dsat(2.5)))
print("Id_sat: " + str(sah_I_D_sat_nmos(V_dsat(2.5)+0.001, 2.5)))

print("Vdsat ratio: " + str(nsat_Vdsat(2.5)/psat_Vdsat(2.5)))
print("Isat ratio: " + str(vsat_sah_I_D_sat_nmos(1, 2.5)/vsat_sah_I_D_sat_pmos(1, 2.5)))