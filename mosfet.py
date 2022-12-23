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
X_ox = 50*np.power(10.0, -8) #A->cm
Na = 4*np.power(10.0, 17)
PHI_m = 4.10
u_nch = 250
W = 15*np.power(10.0,-4) #um cm
L = 5*np.power(10.0,-4)
Vtn = 0.4338
lamb = 0 # Channel length modulation factor

### USEFUL PARAMETERS
C_ox = e_ox/X_ox
n_b_o = n_i_2/Na
p_b_o = Na
Ef_Ev = lambda Na : -kbT*np.log(N_v/Na)
PHI_pSi = lambda Na : X_si + Eg + Ef_Ev(Na)
phi_pm = lambda Na : PHI_m - PHI_pSi(Na)
L_db = lambda Na : np.power((e_s*kbT)/(q*Na), 0.5)
phi_fb = lambda Na : kbT*np.log(Na/n_i)
Q_fit = q*(3*np.power(10.0, 10))

### Flatband, Threshold Voltages
V_FB = lambda Na : phi_pm(Na) - (Q_fit/C_ox)
V_TN = lambda Na : V_FB(Na) + (2*phi_fb(Na)) + (np.sqrt(4*q*e_s*Na*phi_fb(Na))/C_ox)

print("Threshold is: " + p.sn(V_TN(Na)) + " V")

# Mobility Degradation
d_u_nch = lambda Na, theta_n, Vgs : (0.5*mobi.u_n_P(Na))/(1+np.power((((Vgs-Vtn)/(theta_n*4*np.power(10.0, -7)))),1.85))

### Level 1 Sah
Id_lin_1 = lambda Vds, Vgs : u_nch*C_ox*(W/L)*((Vgs-Vtn)*Vds-(0.5*Vds*Vds))*((lamb*Vds)+1) if (Vgs-Vtn)>0 else (Vds-Vds)
Id_sat_1 = lambda Vds, Vgs : u_nch*C_ox*(W/(2*L))*np.power(Vgs-Vtn,2)*(Vds/Vds)*((lamb*Vds)+1) if (Vgs-Vtn)>0 else Vds-Vds
Vd_sat_1 = lambda Vgs : Vgs-Vtn if (Vgs-Vtn)>0 else 0



### Level 2 Ihantola-Moll
Id_lin_2 = lambda Vds, Vgs, Na : u_nch*C_ox*(W/L)*((Vgs-V_FB(Na)-2*phi_fb(Na)-(0.5*Vds))*Vds-(2/3)*((1/C_ox)*(np.sqrt(2*q*e_s*Na))*(np.power(2*phi_fb(Na)+Vds,1.5)-np.power(2*phi_fb(Na),1.5))))*((lamb*Vds)+1)
Vd_sat_2 = lambda Vgs, Na : Vgs - V_FB(Na) - 2*phi_fb(Na) +(q*e_s*Na/(C_ox*C_ox))*(1-np.sqrt(1+((2*C_ox*C_ox)/(q*e_s*Na))*(Vgs-V_FB(Na)),))
Id_sat_2 = lambda Vds, Vgs, Na : u_nch*C_ox*(W/L)*((Vgs-V_FB(Na)-2*phi_fb(Na)-(0.5*Vd_sat_2(Vgs,Na)))*Vd_sat_2(Vgs,Na)-(2/3)*((1/C_ox)*(np.sqrt(2*q*e_s*Na))*(np.power(2*phi_fb(Na)+Vd_sat_2(Vgs,Na),1.5)-np.power(2*phi_fb(Na),1.5))))*((lamb*Vds)+1)

### Level 3 MBC (With mobility degradation)
a = lambda Na : 1+(1/C_ox)*np.power((q*e_s*Na)/(4*phi_fb(Na)), 0.5)
Id_lin_3 = lambda Vds, Vgs : d_u_nch(Na, 3.6*np.power(10.0,6),Vgs)*C_ox*(W/L)*((Vgs-Vtn)*Vds-(a(Na)*0.5*Vds*Vds))*((lamb*Vds)+1)
Id_sat_3 = lambda Vds, Vgs : d_u_nch(Na, 3.6*np.power(10.0,6),Vgs)*C_ox*(W/(2*L))*(1/a(Na))*np.power(Vgs-Vtn,2)*((lamb*Vds)+1)
Vd_sat_3 = lambda Vgs : (1/a(Na))*(Vgs-Vtn)


print("Vdsat: " + p.sn(Vd_sat_1(1.5)))
print("Id at Vds = 1.5 using SAH: " + p.sn(Id_sat_1(1.5, 1.5)))
print("Vdsat using IM: " + p.sn(Vd_sat_2(1.5, Na)))
print("Id at Vds = 1.5 using IM: " + p.sn(Id_sat_2(1.5, 1.5, Na)))
print("Id ration between Sah and IM in sat " + p.sn((Id_sat_1(1.5, 1.5)-Id_sat_2(1.5,1.5, Na))/(Id_sat_2(1.5, 1.5, Na))))

Vds_1_lin_0 = np.linspace(0, Vd_sat_1(0), 1000)
Vds_1_lin_1 = np.linspace(0, Vd_sat_1(1), 1000)
Vds_1_lin_2 = np.linspace(0, Vd_sat_1(2), 1000)
Vds_1_lin_3 = np.linspace(0, Vd_sat_1(3), 1000)
Vds_1_lin_4 = np.linspace(0, Vd_sat_1(4), 1000)

Vds_2_lin_0 = np.linspace(0, Vd_sat_2(0, Na), 1000)
Vds_2_lin_1 = np.linspace(0, Vd_sat_2(1, Na), 1000)
Vds_2_lin_2 = np.linspace(0, Vd_sat_2(2, Na), 1000)
Vds_2_lin_3 = np.linspace(0, Vd_sat_2(3, Na), 1000)
Vds_2_lin_4 = np.linspace(0, Vd_sat_2(4, Na), 1000)

Vds_3_lin_0 = np.linspace(0, Vd_sat_3(0), 1000)
Vds_3_lin_1 = np.linspace(0, Vd_sat_3(1), 1000)
Vds_3_lin_2 = np.linspace(0, Vd_sat_3(2), 1000)
Vds_3_lin_3 = np.linspace(0, Vd_sat_3(3), 1000)
Vds_3_lin_4 = np.linspace(0, Vd_sat_3(4), 1000)

Vds_1_sat_0 = np.linspace(Vd_sat_1(0), 3, 1000)
Vds_1_sat_1 = np.linspace(Vd_sat_1(1), 3, 1000)
Vds_1_sat_2 = np.linspace(Vd_sat_1(2), 3, 1000)
Vds_1_sat_3 = np.linspace(Vd_sat_1(3), 3, 1000)
Vds_1_sat_4 = np.linspace(Vd_sat_1(4), 3, 1000)

Vds_2_sat_0 = np.linspace(Vd_sat_2(0, Na), 1000)
Vds_2_sat_1 = np.linspace(Vd_sat_2(1, Na), 1000)
Vds_2_sat_2 = np.linspace(Vd_sat_2(2, Na), 1000)
Vds_2_sat_3 = np.linspace(Vd_sat_2(3, Na), 1000)
Vds_2_sat_4 = np.linspace(Vd_sat_2(4, Na), 1000)

Vds_3_sat_0 = np.linspace(Vd_sat_3(0), 3, 1000)
Vds_3_sat_1 = np.linspace(Vd_sat_3(1), 3, 1000)
Vds_3_sat_2 = np.linspace(Vd_sat_3(2), 3, 1000)
Vds_3_sat_3 = np.linspace(Vd_sat_3(3), 3, 1000)
Vds_3_sat_4 = np.linspace(Vd_sat_3(4), 3, 1000)


figure, axis = plt.subplots(2, 2)

axis[0,0].plot(Vds_1_lin_0, Id_lin_1(Vds_1_lin_0, 0),)
axis[0,0].plot(Vds_1_lin_1, Id_lin_1(Vds_1_lin_1, 1),'b')
axis[0,0].plot(Vds_1_lin_2, Id_lin_1(Vds_1_lin_2, 2),'g')
axis[0,0].plot(Vds_1_lin_3, Id_lin_1(Vds_1_lin_3, 3),'y')
axis[0,0].plot(Vds_1_lin_4, Id_lin_1(Vds_1_lin_4, 4),'r')
axis[0,0].plot(Vds_1_sat_0, Id_sat_1(Vds_1_sat_0, 0))
axis[0,0].plot(Vds_1_sat_1, Id_sat_1(Vds_1_sat_1, 1),'b')
axis[0,0].plot(Vds_1_sat_2, Id_sat_1(Vds_1_sat_2, 2),'g')
axis[0,0].plot(Vds_1_sat_3, Id_sat_1(Vds_1_sat_3, 3),'y')
axis[0,0].plot(Vds_1_sat_4, Id_sat_1(Vds_1_sat_4, 4),'r')
axis[0,0].set_xlabel("$V_{DS}$ V")
axis[0,0].set_ylabel("$I_D(V_{DS})$ A")
axis[0,0].grid()
axis[0,0].legend(['$V_{gs}$=0', '$V_{gs}$=1', '$V_{gs}$=2', '$V_{gs}$=3', '$V_{gs}$=4'])
axis[0,0].set_title('$I_{D}$ Level 1 - Sah Model')

axis[0,1].plot(Vds_3_lin_0, Id_lin_3(Vds_3_lin_0, 0.001))
axis[0,1].plot(Vds_3_lin_1, Id_lin_3(Vds_3_lin_1, 1),'b')
axis[0,1].plot(Vds_3_lin_2, Id_lin_3(Vds_3_lin_2, 2),'g')
axis[0,1].plot(Vds_3_lin_3, Id_lin_3(Vds_3_lin_3, 3),'y')
axis[0,1].plot(Vds_3_lin_4, Id_lin_3(Vds_3_lin_4, 4),'r')
axis[0,1].plot(Vds_3_sat_0, Id_sat_3(Vds_3_sat_0, 0.001))
axis[0,1].plot(Vds_3_sat_1, Id_sat_3(Vds_3_sat_1, 1),'b')
axis[0,1].plot(Vds_3_sat_2, Id_sat_3(Vds_3_sat_2, 2),'g')
axis[0,1].plot(Vds_3_sat_3, Id_sat_3(Vds_3_sat_3, 3),'y')
axis[0,1].plot(Vds_3_sat_4, Id_sat_3(Vds_3_sat_4, 4),'r')
axis[0,1].set_xlabel("$V_{DS}$ V")
axis[0,1].set_ylabel("$I_D(V_{DS})$ A")
axis[0,1].grid()
axis[0,1].legend(['$V_{gs}$=0', '$V_{gs}$=1', '$V_{gs}$=2', '$V_{gs}$=3', '$V_{gs}$=4'])
axis[0,1].set_title('Static Drain Current Characteristics Lvelel 3')


figure.tight_layout();


#plt.show()
