# importing the modules
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

### Problem parameters
X_ox = 50*np.power(10.0, -8) #A->cm (1 nm = 10 A)
Na = 4*np.power(10.0, 17)
PHI_m = 4.10
u_nch = 250
W = 15*np.power(10.0,-4) #um
L = 0.065*np.power(10.0,-4) #um
Vtn = 0.2181
lamb = 0 # Channel length modulation factor
e_sat = 3.8*np.power(10.0, 4)
v_sat = u_nch*0.5*e_sat

### USEFUL PARAMETERS
C_ox = e_ox/X_ox
n_b_o = n_i_2/Na
p_b_o = Na
Ef_Ev = lambda Na : -kbT*np.log(N_v/Na)
PHI_pSi = lambda Na : X_si + Eg + Ef_Ev(Na)
phi_pm = lambda Na : PHI_m - PHI_pSi(Na)
L_db = lambda Na : np.power((e_s*kbT)/(q*Na), 0.5)
phi_fb = lambda Na : kbT*np.log(Na/n_i)
Q_fit = q*(2*np.power(10.0, 10))

### Flatband, Threshold Voltages
V_FB = lambda Na : phi_pm(Na) - (Q_fit/C_ox)
V_TN = lambda Na : V_FB(Na) + (2*phi_fb(Na)) + (np.sqrt(4*q*e_s*Na*phi_fb(Na))/C_ox)
print("Threshold is: " + p.sn(V_TN(Na)) + " V")
print("SC Threshold is: " + p.sn(V_TN(Na)/2))

### Level 1 Sah - YES VSAT
Id_lin_1 = lambda Vds, Vgs : u_nch*C_ox*(W/L)*(1/(1+(Vds/(e_sat*L))))*((Vgs-Vtn)*Vds-0.5*Vds*Vds) if Vgs>Vtn else 0*Vds
Id_sat_1 = lambda Vds, Vgs : v_sat*C_ox*W*(np.power(Vgs-Vtn, 2)/((Vgs-Vtn)+e_sat*L))*((Vds+1)/(Vds+1)) if Vgs>Vtn else 0*Vds
Vd_sat_1 = lambda Vgs : ((Vgs-Vtn)*e_sat*L)/((Vgs-Vtn)+L*e_sat) if Vgs > 0 else 0*(Vgs)
print("Vdsat: " + p.sn(Vd_sat_1(1.5)))
print("Id(Vdsat): " + p.sn(Id_sat_1(1.5, 1.5)))


### Level 1 Sah - NO VSAT
Id_lin_2 = lambda Vds, Vgs : u_nch*C_ox*(W/L)*((Vgs-Vtn)*Vds-(0.5*Vds*Vds))*((lamb*Vds)+1) if (Vgs-Vtn)>0 else (Vds-Vds)
Id_sat_2 = lambda Vds, Vgs : u_nch*C_ox*(W/(2*L))*np.power(Vgs-Vtn,2)*(Vds/Vds)*((lamb*Vds)+1) if (Vgs-Vtn)>0 else Vds-Vds
Vd_sat_2 = lambda Vgs : Vgs-Vtn if (Vgs-Vtn)>0 else 0
print("Vdsat noVSAT: " + p.sn(Vd_sat_2(1.5)))
print("Id(Vdsat) noVSAT: " + p.sn(Id_sat_2(1.5, 1.5)))

Vds_1_lin_0 = np.linspace(0, Vd_sat_1(0), 1000)
Vds_1_lin_1 = np.linspace(0, Vd_sat_1(1), 1000)
Vds_1_lin_2 = np.linspace(0, Vd_sat_1(2), 1000)
Vds_1_lin_3 = np.linspace(0, Vd_sat_1(3), 1000)
Vds_1_lin_4 = np.linspace(0, Vd_sat_1(4), 1000)
Vds_1_sat_0 = np.linspace(Vd_sat_1(0), 5, 1000)
Vds_1_sat_1 = np.linspace(Vd_sat_1(1), 5, 1000)
Vds_1_sat_2 = np.linspace(Vd_sat_1(2), 5, 1000)
Vds_1_sat_3 = np.linspace(Vd_sat_1(3), 5, 1000)
Vds_1_sat_4 = np.linspace(Vd_sat_1(4), 5, 1000)

Vds_2_lin_0 = np.linspace(0, Vd_sat_2(0), 1000)
Vds_2_lin_1 = np.linspace(0, Vd_sat_2(1), 1000)
Vds_2_lin_2 = np.linspace(0, Vd_sat_2(2), 1000)
Vds_2_lin_3 = np.linspace(0, Vd_sat_2(3), 1000)
Vds_2_lin_4 = np.linspace(0, Vd_sat_2(4), 1000)
Vds_2_sat_0 = np.linspace(Vd_sat_2(0), 5, 1000)
Vds_2_sat_1 = np.linspace(Vd_sat_2(1), 5, 1000)
Vds_2_sat_2 = np.linspace(Vd_sat_2(2), 5, 1000)
Vds_2_sat_3 = np.linspace(Vd_sat_2(3), 5, 1000)
Vds_2_sat_4 = np.linspace(Vd_sat_2(4), 5, 1000)

figure, axis = plt.subplots(2)
axis[0].plot(Vds_1_lin_0, Id_lin_1(Vds_1_lin_0, 0))
axis[0].plot(Vds_1_lin_1, Id_lin_1(Vds_1_lin_1, 1), 'r')
axis[0].plot(Vds_1_lin_2, Id_lin_1(Vds_1_lin_2, 2),'y')
axis[0].plot(Vds_1_lin_3, Id_lin_1(Vds_1_lin_3, 3),'g')
axis[0].plot(Vds_1_lin_4, Id_lin_1(Vds_1_lin_4, 4),'b')
axis[0].plot(Vds_1_sat_0, Id_sat_1(Vds_1_sat_0, 0))
axis[0].plot(Vds_1_sat_1, Id_sat_1(Vds_1_sat_1, 1),'r')
axis[0].plot(Vds_1_sat_2, Id_sat_1(Vds_1_sat_2, 2),'y')
axis[0].plot(Vds_1_sat_3, Id_sat_1(Vds_1_sat_3, 3),'g')
axis[0].plot(Vds_1_sat_4, Id_sat_1(Vds_1_sat_4, 4),'b')
axis[0].set_xlabel("$V_{DS}$ V")
axis[0].set_ylabel("$I_D(V_{DS})$ A")
axis[0].grid()
axis[0].legend(['$V_{gs}$=0', '$V_{gs}$=1', '$V_{gs}$=2', '$V_{gs}$=3', '$V_{gs}$=4'])
axis[0].set_title('Short-Channel $I_{D}$ Characteristics Level 1 with VSat')
axis[1].plot(Vds_2_lin_0, Id_lin_2(Vds_2_lin_0, 0))
axis[1].plot(Vds_2_lin_1, Id_lin_2(Vds_2_lin_1, 1), 'r')
axis[1].plot(Vds_2_lin_2, Id_lin_2(Vds_2_lin_2, 2),'y')
axis[1].plot(Vds_2_lin_3, Id_lin_2(Vds_2_lin_3, 3),'g')
axis[1].plot(Vds_2_lin_4, Id_lin_2(Vds_2_lin_4, 4),'b')
axis[1].plot(Vds_2_sat_0, Id_sat_2(Vds_2_sat_0, 0))
axis[1].plot(Vds_2_sat_1, Id_sat_2(Vds_2_sat_1, 1),'r')
axis[1].plot(Vds_2_sat_2, Id_sat_2(Vds_2_sat_2, 2),'y')
axis[1].plot(Vds_2_sat_3, Id_sat_2(Vds_2_sat_3, 3),'g')
axis[1].plot(Vds_2_sat_4, Id_sat_2(Vds_2_sat_4, 4),'b')
axis[1].set_xlabel("$V_{DS}$ V")
axis[1].set_ylabel("$I_D(V_{DS})$ A")
axis[1].grid()
axis[1].legend(['$V_{gs}$=0', '$V_{gs}$=1', '$V_{gs}$=2', '$V_{gs}$=3', '$V_{gs}$=4'])
axis[1].set_title('Short-Channel $I_{D}$ Characteristics Level 1 without VSat')

plt.show()