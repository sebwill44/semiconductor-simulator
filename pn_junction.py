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

### Parameters
Na = 9*np.power(10.0, 19)
Nd = 2.5*np.power(10.0, 15)
Vpn = 0
Wp = 0
Nd2 = 7.5*np.power(10.0, 16)


V_bi = lambda Na, Nd, Vpn : kbT*np.log((Na*Nd)/n_i_2)-Vpn
print("V_bi: " + p.sn(V_bi(Na,Nd,Vpn)) + " V")

### Charges, E-Field, Voltage
sigma_P = lambda Na : -q*Na
sigma_N = lambda Nd : q*Nd

E_P = lambda Na, Nd, x, Vpn : (-q*Na/e_s)*(x+W_dP(Na,Nd,Vpn))
E_N = lambda Na, Nd, x, Vpn : (-q*Nd/e_s)*(W_dN(Na,Nd,Vpn)-x)

vbi_P = lambda Na, Nd, Vpn : (Nd/(Na+Nd))*(V_bi(Na,Nd,Vpn))
vbi_N = lambda Na, Nd, Vpn : (Na/(Na+Nd))*(V_bi(Na,Nd,Vpn))

phi_P = lambda Na, Nd, x, Vpn : (q*Na/(2*e_s))*np.power((x+W_dP(Na,Nd,Vpn)),2)
phi_N = lambda Na, Nd, x, Vpn : V_bi(Na,Nd,Vpn)-(q*Nd/(2*e_s))*np.power((W_dN(Na,Nd,Vpn)-x),2)

W_dP = lambda Na, Nd, Vpn : np.sqrt((2*e_s/q)*Nd*(V_bi(Na,Nd,Vpn))/(Na*(Na+Nd)))
W_dN = lambda Na, Nd, Vpn : np.sqrt((2*e_s/q)*Na*(V_bi(Na,Nd,Vpn))/(Nd*(Na+Nd)))
W_d = lambda Na, Nd, Vpn : W_dP(Na,Nd,Vpn) + W_dN(Na,Nd,Vpn)
print("Depletion width total: " + p.sn(W_d(Na, Nd, Vpn)) + " cm = " + p.sn(W_d(Na,Nd,Vpn)*np.power(10.0, 4)) + " um")

### Calculate total length
Wn = 0

### P SIDE
n_P_o = lambda Na, Nd, x, Vpn : np.exp(Vpn/kbT)*(n_i_2/Na)*np.exp((q*Na/(2*e_s*kbT))*np.power((W_dP(Na,Nd,Vpn)+x),2))
p_P_o = lambda Na, Nd, x, Vpn : Na*np.exp((-q*Na/(2*e_s*kbT))*np.power((W_dP(Na,Nd,Vpn)+x),2))

### N SIDE
n_N_o = lambda Na, Nd, x, Vpn : Nd*np.exp((-q*Nd/(2*e_s*kbT))*np.power((W_dN(Na,Nd,Vpn)-x),2))
p_N_o = lambda Na, Nd, x, Vpn : np.exp(Vpn/kbT)*(n_i_2/Nd)*np.exp((q*Nd/(2*e_s*kbT))*np.power((W_dN(Na,Nd,Vpn)-x),2))

### Long Base N Side quasi neutral region
p_N_q_lb = lambda Na, Nd, x, Vpn : (n_i_2/Nd)*np.exp(-(x-W_dN(Na,Nd,Vpn))/mobi.L_p_N(Nd))*(np.exp(Vpn/kbT)-1)+((n_i_2/Nd))

### Short Base N Side Quasi Neutral Region
p_N_q_sb = lambda Na, Nd, x, Vpn : (n_i_2/Nd)*np.exp((Wn-x)/(Wn-W_dN(Na,Nd,Vpn)))*(np.exp(Vpn/kbT)-1)+((n_i_2/Nd))

### Long Base P Side Quasi Neutral Region
n_P_q_lb = lambda Na, Nd, x, Vpn : (n_i_2/Na)*np.exp((x+W_dP(Na,Nd,Vpn))/mobi.L_n_P(Na))*(np.exp(Vpn/kbT)-1)+((n_i_2/Na))

### Short Base P Side Quasi Neutral Region
n_P_q_sb = lambda Na, Nd, x, Vpn : (n_i_2/Na)*np.exp((Wp+x)/(Wp-W_dP(Na,Nd,Vpn)))*(np.exp(Vpn/kbT)-1)+((n_i_2/Na))



### Current Equations
J_p_diff_N_lb = lambda Na, Nd, x, Vpn : ((q*mobi.D_p_N(Nd)*n_i_2)/(mobi.L_p_N(Nd)*Nd))*np.exp(-(x-W_dN(Na,Nd,Vpn))/(mobi.L_p_N(Nd)))*(np.exp(Vpn/kbT)-1)
J_n_diff_P_lb = lambda Na, Nd, x, Vpn : ((q*mobi.D_n_P(Na)*n_i_2)/(mobi.L_n_P(Na)*Na))*np.exp((x+W_dP(Na,Nd,Vpn))/(mobi.L_n_P(Na)))*(np.exp(Vpn/kbT)-1)

ratio = lambda Na, Nd, Nd2 : ((mobi.D_n_P(Na)/(mobi.L_n_P(Na)*Na))+(mobi.D_p_N(Nd)/(mobi.L_p_N(Nd)*Nd)))*(1/((mobi.D_n_P(Na)/(mobi.L_n_P(Na)*Na))+(mobi.D_p_N(Nd2)/(mobi.L_p_N(Nd2)*Nd2))))
print("Ration of diff currents : " + p.sn(ratio(Na, Nd, Nd2)))

J_p_diff_N_sb = lambda Na, Nd, x, Vpn : ((q*mobi.D_p_N(Nd)*n_i_2)/((Wn-W_dN(Na,Nd,Vpn))*Nd))*(np.exp(Vpn/kbT)-1)
J_n_diff_P_sb = lambda Na, Nd, x, Vpn : ((q*mobi.D_n_P(Na)*n_i_2)/(mobi.L_n_P(Na)*Na))*(np.exp(Vpn/kbT)-1)

J_D_scr = lambda Na, Nd, Vpn : ((q*n_i*W_d(Na,Nd,Vpn))/(recombo.T_rec_P(Na)+recombo.T_rec_N(Nd)))*(np.exp(Vpn/(2*kbT))-1)
print("Ration of SCR currents : " + p.sn(J_D_scr(Na, Nd, 0.5)/J_D_scr(Na, Nd2, 0.5)))



J_D = lambda Na, Nd, Vpn : J_n_diff_P_lb(Na, Nd, W_dP(Na, Nd, Vpn),Vpn) + J_p_diff_N_lb(Na, Nd, W_dN(Na,Nd,Vpn),Vpn) + J_D_scr(Na, Nd, Vpn)
### BONUS Q
print("this should be 20: " + p.sn(J_D_scr(Na, Nd, .33)/(J_D(Na, Nd,.33))))
### BONUS Q

J_n_N = lambda Na, Nd, Vpn : J_D(Na, Nd, Vpn) - J_p_diff_N_lb(Na, Nd, W_dP(Na,Nd,Vpn),Vpn)
J_p_P = lambda Na, Nd, Vpn : J_D(Na, Nd, Vpn) - J_n_diff_P_lb(Na,Nd, W_dN(Na,Nd,Vpn), Vpn)

print("Hole diffusion current at nside sb: " + p.sn(J_p_diff_N_sb(Na,Nd,W_dN(Na,Nd,Vpn),Vpn)*(2*np.power(10.0,-4))) + " A/cm^2")
print("Electron diffusion current at pside sb: " + p.sn(J_n_diff_P_sb(Na,Nd,W_dP(Na,Nd,Vpn),Vpn)*(2*np.power(10.0,-4))) + " A/cm^2")
print("Space charge recombo current: " + p.sn(J_D_scr(Na,Nd,Vpn)*(2*np.power(10.0,-4))) + " A/cm^2")
print(p.sn(J_n_N(Na,Nd,Vpn)))

### Capacitances (long base)
C_pn_dep = lambda Na, Nd, Vpn : np.sqrt((q*e_s/2)*((Na*Nd)/(Na + Nd))*(1/(V_bi(Na,Nd,Vpn))))
C_pn_N_diff_lb = lambda Na, Nd : ((q*n_i_2*mobi.L_p_N(Nd))/(kbT*Nd))
C_pn_P_diff_lb = lambda Na, Nd : ((q*n_i_2*mobi.L_n_P(Na))/(kbT*Na))
C_pn_diff = lambda Na, Nd, Vpn : (C_pn_N_diff_lb(Na, Nd) + C_pn_P_diff_lb(Na, Nd))*np.exp(Vpn/kbT)
C_pn = lambda Na, Nd, Vpn : C_pn_diff(Na,Nd,Vpn) + C_pn_dep(Na,Nd,Vpn)

print("Dep cap ration :" + p.sn(C_pn_dep(Na, Nd, 0)/C_pn_dep(Na, Nd2, 0)))

print(p.sn(C_pn_dep(Na,Nd,Vpn)))
print(p.sn(C_pn_diff(Na,Nd,Vpn)))

px = np.linspace(-W_dP(Na,Nd,Vpn), 0, 1000)
nx = np.linspace(0, W_dN(Na,Nd,Vpn), 1000)
qnx = np.linspace(W_dN(Na,Nd,Vpn), W_dN(Na,Nd,Vpn)+Wn, 1000)
qpx = np.linspace(-(Wp+W_dP(Na,Nd,Vpn)), -W_dP(Na,Nd,Vpn), 1000)

print((recombo.T_rec_P(Na)))
fig, axis = plt.subplots(2,2)

axis[1,1].plot(px, E_P(Na,Nd,px,Vpn), 'r')
axis[1,1].plot(nx, E_N(Na,Nd,nx,Vpn), 'r')
axis[1,0].plot(px, phi_P(Na,Nd,px,Vpn), 'r')
axis[1,0].plot(nx, phi_N(Na,Nd,nx,Vpn),'r')
axis[1,0].set_title("Potential in Depletion Region")
axis[1,0].set_ylabel(r"Potential, $V$")
axis[1,0].set_xlabel("Position $x$ ($cm$) ")
axis[1,0].grid()

axis[1,1].set_title("E-field in Depletion Region")
axis[1,1].set_ylabel(r"E-field, $\frac{V}{cm}$")
axis[1,1].set_xlabel("Position $x$ ($cm$) ")
axis[1,1].grid()
### Uncomment below for E field, V, and charge distributions
""" axis[0].plot(qnx, p_N_q_lb(Na,Nd,qnx,Vpn),'r')
axis[0].plot(qpx, n_P_q_lb(Na,Nd,nx,Vpn),'g') """

""" axis[0].plot(px, E_P(Na,Nd,px,Vpn))
axis[0].plot(nx, E_N(Na,Nd,nx,Vpn))

axis[1].plot(px, phi_P(Na,Nd,px,Vpn))
axis[1].plot(nx, phi_N(Na,Nd,nx,Vpn))

axis[2].plot(px, sigma_P(Na)*(px/px))
axis[2].plot(nx, sigma_N(Nd)*(nx/nx)) """


axis[0,0].plot(px, n_P_o(Na,Nd,px,Vpn), 'g')
axis[0,0].plot(px, p_P_o(Na,Nd,px,Vpn),'r')
axis[0,0].plot(nx, n_N_o(Na,Nd,nx,Vpn),'g')
axis[0,0].plot(nx, p_N_o(Na,Nd,nx,Vpn),'r')
axis[0,0].set_yscale('log')
axis[0,0].set_title("Electron- and Hole-Concentration Distributions")
axis[0,0].set_ylabel("Carrier-Concentration Distributions, $n^{o}(x)$ and $p^{o}(x)$ ($cm^{-3}$)")
axis[0,0].set_xlabel("Position $x$ ($cm$) ")
axis[0,0].grid()
axis[0,0].legend(['Electrons', 'Holes'])

bv = np.linspace(0.1, 0.5, 1000)
axis[0,1].plot(bv, C_pn_dep(Na,Nd,bv), 'g')
axis[0,1].plot(bv, C_pn_diff(Na,Nd,bv),'r')
axis[0,1].plot(bv, C_pn(Na,Nd,bv),'y')
axis[0,1].set_title("PN Junction Capacitance vs. Bias Voltage")
axis[0,1].set_ylabel("Capacitance pF")
axis[0,1].set_xlabel(r"$V_{PN} (V)$")
axis[0,1].legend(['Depletion', 'Diffusion', 'Total'])
axis[0,1].grid()

fig.tight_layout()
#plt.show()