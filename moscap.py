import numpy as np
import matplotlib.pyplot as plt
import printing as p
import mobilities as mobi
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

### PROBLEM PARAMETERS
Na = 4*np.power(10.0, 17)
PHI_m = 4.10
X_ox = 50*np.power(10.0, -8) #A->cm
Q_fit = q*(3*np.power(10.0, 10))
Q_scb_o = 1.51008*np.power(10.0,-7)
V_scb_o = 0.71254

### USEFUL PARAMETERS
C_ox = e_ox/X_ox
n_b_o = n_i_2/Na
p_b_o = Na
Ef_Ev = lambda Na : -kbT*np.log(N_v/Na)
PHI_pSi = lambda Na : X_si + Eg + Ef_Ev(Na)
phi_pm = lambda Na : PHI_m - PHI_pSi(Na)
L_db = lambda Na : np.power((e_s*kbT)/(q*Na), 0.5)
phi_fb = lambda Na : kbT*np.log(Na/n_i)

### Flatband, Threshold Voltages
V_FB = lambda Na : phi_pm(Na) - (Q_fit/C_ox)
V_TN = lambda Na : V_FB(Na) + (2*phi_fb(Na)) + (np.sqrt(4*q*e_s*Na*phi_fb(Na))/C_ox)

### Tes problem
print("2 Phi fb: " +p.sn(2*phi_fb(Na)))
print("PHI_pm for this problem: " + p.sn(phi_pm(Na)) + " eV")
print("Vox(Vtn): " + p.sn(2*phi_fb(Na) + phi_pm(Na) - V_TN(Na)))
# Oxide voltage, field
E_ox = -(Q_scb_o-Q_fit)/e_ox
V_ox = X_ox*E_ox
E_ox_vtn = lambda Na : (1/e_ox)*(np.sqrt(4*q*e_s*Na*phi_fb(Na))-(Q_fit))
V_ox_vtn = lambda Na : X_ox*(E_ox_vtn(Na))
V_ox_vgb = lambda Na : 5*V_ox_vtn(Na)
E_ox_vgb = lambda Na : V_ox_vgb(Na)/(X_ox)
print("E_ox_vgb " + p.sn(E_ox_vgb(Na))) 


### Printing
print("Vfb: " + p.sn(V_FB(Na)) + " V")
print("Vtn: " + p.sn(V_TN(Na)) + " V")
print("Cox = " + p.sn(C_ox) + " F/cm^2")
print("PHI_pSi = " + p.sn(PHI_pSi(Na)) + " eV")
print("PHI_pm = " + p.sn(phi_pm(Na)) + " eV")
print("MAKE SURE TO PLUG IN CORRECT QSCBE_ox: " + p.sn(E_ox) + "V/cm")
print("Vox: " + p.sn(V_ox) + " V")

### Key equations
def Q_scb_1(V_scb):
    if (C_ox*phi_pm(Na) - Q_fit) < 0:
        print("Q_scb0 is negative! Charge is NEGATIVE!")
        return -1*(C_ox*V_scb + C_ox*phi_pm(Na) - Q_fit)
    else:
        return C_ox*V_scb + C_ox*phi_pm - Q_fit

def Q_scb_2(V_scb):
    return (root_2*e_s*kbT/(L_db(Na)))*np.power(((np.exp(-B*V_scb)+B*V_scb-1)+(n_b_o/p_b_o)*(np.exp(B*V_scb)-B*V_scb-1)),0.5)

print("E" + p.sn(((E_ox_vgb(Na)*e_s)+(Q_fit))/(e_s)))
print("Q" + p.sn(e_s*(((E_ox_vgb(Na)*e_s)+(Q_fit))/(e_s))))
print("Vgb :" + p.sn(phi_pm(Na)+1.17+((1/C_ox)*Q_fit+5.122*np.power(10.0, -6))))
### Stuff based off Vscb and Qscb
V_scb = np.linspace(-2, 2, 1000)
plt.plot(V_scb, Q_scb_2(V_scb))
plt.plot(V_scb, Q_scb_1(V_scb))
#axis[0,0].set_xscale('log')
plt.yscale('log')
plt.title("Substrate Space Charge Density vs. Substrate Space Charge Voltage")
plt.xlabel("Substrate Space Charge Voltage, $V_{SC,B}$(V)")
plt.ylabel("Substrate Space Charge Density $|Q_{SC,B}(V_{SC,B})|$ (C/cm$^{2}$)")
plt.grid()

plt.show()
