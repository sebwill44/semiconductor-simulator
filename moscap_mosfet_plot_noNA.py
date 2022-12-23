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
theta_n = 3.6*np.power(10.0, 6)

def phi_si(Na):
    return X_si + Eg + kbT*np.log(Na/N_v)

def phi_pm(Na):
    return 4.125 - phi_si(Na)

def phi_fb(Na):
    return kbT*np.log(Na/n_i)

t_ox = 4*np.power(10.0, -7)
C_ox = epsilon_sio2/t_ox
print(str(C_ox))
V_ox_fb = -(2*np.power(10.0, 10)*q)/C_ox
print(V_ox_fb)

def V_fb(Na):
    return V_ox_fb + phi_pm(Na)

a = 1+(1/C_ox)*np.power((q*epsilon_s*2.3*np.power(10.0, 17))/(4*phi_fb(2.3*np.power(10.0, 17))), 0.5)

print("a:" + str(a))
def u_nP(Na):
    return 232 + (1180)/(1+np.power(Na/(8*np.power(10.0, 16)), 0.9))
print(u_nP(2.6*np.power(10.0, 17)))
# Finding NA
def V_tn(Na):
    return -0.25 + V_fb(Na) + 2*phi_fb(Na) + (1/C_ox)*np.power(4*q*epsilon_s*Na*phi_fb(Na),0.5)

def u_nCh(Vgs):
    return (0.5*u_nP(2.6*np.power(10.0,17)))/(1+np.power((((Vgs-0.25)/(theta_n*4*np.power(10.0, -7)))),1.85))

print(str(u_nCh(1)))
print(str(u_nCh(2)))
print(str(u_nCh(3)))
print(str(u_nCh(4)))

Vdsat_1 = (1-0.25)/a
Vdsat_2 = (2-0.25)/a
Vdsat_3 = (3-0.25)/a
Vdsat_4 = (4-0.25)/a

def I_D_lin(Vds, Vgs):
    return u_nCh(Vgs)*C_ox*10*(((Vgs-0.25)*Vds)-((0.5*a*Vds*Vds)))*(1+0.04*(Vds-((Vgs-0.25)/a)))

def I_D_sat(Vds, Vgs):
    return u_nCh(Vgs)*C_ox*5*np.power(Vgs-0.25, 2)*(1/(a))*(Vds/Vds)*(1+0.04*(Vds-((Vgs-0.25)/a)))

def I_D_sth(Vds, Vgs):
    return u_nCh(Vgs)*C_ox*10*(a-1)*(kbT*kbT)*np.exp(q*(Vgs-0.25)/(a*kbT))*(1-np.exp(-Vds/kbT))

B_n = u_nP(2.6*np.power(10.0, 17))*C_ox*0.001
def Vco(Vds, Vgs, y):
    return ((Vgs-0.25)-np.power(np.power(Vgs-0.25, 2)-(1/B_n)*(2*a*I_D_lin(Vds, Vgs)*y),0.5))/(a)

Vds = np.linspace(0, 2, 1000)
Na = np.linspace(np.power(10.0, 12), np.power(10.0, 19))
print("u_nP: " + str(u_nP(2.36*np.power(10.0,17))))

figure, axis = plt.subplots(2, 2)

#axis[0,0].set_xscale('log')
axis[0,0].plot(Na, V_tn(Na), 'r')
axis[0,0].set_xlabel("$N_{a}^{-}$")
axis[0,0].set_ylabel("$V_{TN}$-$V_{TN}(N_{a}^{-})$")
axis[0,0].grid()



Vdss = np.linspace(0,5, 1000)

Vds_1_lin = np.linspace(0, Vdsat_1, 1000)
Vds_2_lin = np.linspace(0, Vdsat_2, 1000)
Vds_3_lin = np.linspace(0, Vdsat_3, 1000)
Vds_4_lin = np.linspace(0, Vdsat_4, 1000)

Vds_1_sat = np.linspace(Vdsat_1, 5, 1000)
Vds_2_sat = np.linspace(Vdsat_2, 5, 1000)
Vds_3_sat = np.linspace(Vdsat_3, 5, 1000)
Vds_4_sat = np.linspace(Vdsat_4, 5, 1000)


#axis[0,1].plot(Vds, I_D_lin(Vds, 0), 'r')
axis[0,1].plot(Vds_1_lin, 1000*I_D_lin(Vds_1_lin, 1), 'g')
axis[0,1].plot(Vds_2_lin, 1000*I_D_lin(Vds_2_lin, 2), 'b')
axis[0,1].plot(Vds_3_lin, 1000*I_D_lin(Vds_3_lin, 3), 'y')
axis[0,1].plot(Vds_4_lin, 1000*I_D_lin(Vds_4_lin, 4), 'r')

#axis[0,1].plot(Vds, I_D_sat(Vds, 0), 'r')
axis[0,1].plot(Vds_1_sat, 1000*I_D_sat(Vds_1_sat, 1), 'g')
axis[0,1].plot(Vds_2_sat, 1000*I_D_sat(Vds_2_sat, 2), 'b')
axis[0,1].plot(Vds_3_sat, 1000*I_D_sat(Vds_3_sat, 3), 'y')
axis[0,1].plot(Vds_4_sat, 1000*I_D_sat(Vds_4_sat, 4), 'r')
axis[0,1].set_xlabel("$V_{DS}$ V")
axis[0,1].set_ylabel("$I_D(V_{DS})$ mA")
axis[0,1].grid()
axis[0,1].legend(['$V_{gs}$=1', '$V_{gs}$=2', '$V_{gs}$=3', '$V_{gs}$=4'])
axis[0,1].set_title('Static Drain Current Characteristics')

y = np.linspace(0, 0.0001)

axis[1,0].plot(y, Vco(1.175, 3, y), 'g')
axis[1,0].set_xlabel("$y$ cm")
axis[1,0].set_ylabel("$V_{CO}(y)$ V")
axis[1,0].grid()

#axis[1,1].set_yscale('log')
figure.tight_layout();


plt.show()