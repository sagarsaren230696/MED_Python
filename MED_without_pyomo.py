import numpy as np
from scipy.optimize import fsolve
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

######################################
# ----------- Effect -----------------
######################################
## Seawater enthalpy parameters 
a = [
    996.7767,
    -3.2406,
    0.0127,
    -4.7723e-5,
    -1.1748,
    0.01169,
    -2.6185e-5,
    7.0661e-8
]
b = [-2.34825e4, 
    3.15183e5, 
    2.80269e6, 
    -1.44606e7, 
    7.82607e3, 
    -4.41733e1, 
    2.1394e-1, 
    -1.99108e4, 
    2.77846e4, 
    9.72801e1
]
c = [
    -4.5838530457E-04,
    2.8230948284E-01,
    1.7945189194E+01,
    1.5361752708E-04,
    5.2669058133E-02,
    6.5604855793E+00,
]

def h_w(T):
    return 141.355 + 4202.07 * (T - 273.15) - 0.535 * (T - 273.15)**2 + 0.004 * (T - 273.15)**3

def h_sw0(T,X):
    hw = h_w(T)
    return hw - X * (
                    -b[0] + b[1]*X + b[2]*X**2 + b[3]*X**3 
                    + b[4] * (T - 273.15) + b[5] * (T - 273.15)**2 + b[6] * (T - 273.15) ** 3
                    + b[7] * X * (T - 273.15) + b[8] * X**2 * (T - 273.15) + b[9] * X * (T - 273.15)**2) 

def h_sw(T,P,X):
    hsw0 = h_sw0(T,X)
    return hsw0 + (P - 0.101) * (a[0] + a[1] * (T - 273.15) + a[2] * (T - 273.15)**2 + a[3] * (T - 273.15)**3 
                                 + X * 1000 * (a[4] + a[5] * (T - 273.15) + a[6] * (T - 273.15)**2 + a[7] * (T - 273.15)**3))

def BPE_sw(T,X):    
    T_degC = T-273.15
    A = c[0] * T_degC**2 + c[1] * T_degC +c[2]
    B = c[3] * T_degC**2 + c[4] * T_degC + c[5]
    BPE = A*X**2+B*X
    return BPE
    
F = 2.0
F1 = F
D = 0.171
X_F = 0.042

T_E = 73.41 + 273.15
T_F = 333.15 #66.81 + 273.15
T_S = 353.15

def Effect_massBalance_1(x):
    B, X_B = x
    eq1 = F - (B + D)
    eq2 = F*X_F - B*X_B
    return [eq1, eq2]

initial_guess = [1.0, 0.05]

sol = fsolve(Effect_massBalance_1,initial_guess)

B, X_B = sol

h_D = CP.PropsSI('H','T',T_E,'Q',1,'Water')
hfg_S = CP.PropsSI('H','T',T_S,'Q',1,'Water') - CP.PropsSI('H','T',T_S,'Q',0,'Water')

P_E = CP.PropsSI('P','T',T_E,'Q',1,'Water') / 1e6
print("Effect pressure is ", P_E)

h_F = h_sw(T_F,P_E,X_F)
h_B = h_sw(T_E,P_E,X_B)

def Effect_energyBalance_1(x):
    m_s = x
    eq = m_s * hfg_S + F * h_F - D * h_D - B * h_B
    return eq

initial_guess_energy = 2.0

sol_energy = fsolve(Effect_energyBalance_1,initial_guess_energy)

m_s, = sol_energy

print("Mass flow rate of steam is ", m_s, " kg/s" )
m_s_1 = m_s 

## Area calculation 
U_E = (1939.1 + 1.40562 * (T_S - 273.15) - 0.02075255 * (T_S - 273.15)**2 + 0.0023186 * (T_S - 273.15)**3)
# print("Heat transfer coefficient: ",U_E)

def Effect_area(x):
    A = x
    eq = m_s*hfg_S - A * U_E * (T_S - T_E)
    return eq

initial_guess_area = 200

sol_area = fsolve(Effect_area, initial_guess_area)

A_E, = sol_area
print("Area of Effect 1 is ", A_E, "m2")

#######################################################
######### Effect 2 onwards ----------------------------
#######################################################
N_effects = 8
net_distillate = m_s_1
T_E_arr = [T_E]
T_S_arr = [T_S]
F_arr = [F]
P_E_arr = [P_E]
m_s_arr = [m_s]

print(f"Effect 1")
print("************")
print("m_s = ",m_s)
print("D = ", D)
print("B = ", B)
print("X_B = ", X_B)
print("T_E = ", T_E-273.15)
print("BPE = ", BPE_sw(T_E, X_B), " K")
print("Net distillate produced = ",net_distillate)
print("")

for N_E in range(2,N_effects+1):
    F = B 
    m_s = D
    T_S = T_E #- BPE_sw(T_E, X_B)
    T_F = T_E
    X_F = X_B

    hfg_S = CP.PropsSI('H','T',T_E,'Q',1,'Water') - CP.PropsSI('H','T',T_S,'Q',0,'Water')

    def Effect_mass_energyBalance(x):
        D, B, X_B, T_E = x
        P_E = CP.PropsSI('P','T',T_E,'Q',1,'Water') / 1e6
        h_F = h_sw(T_F,P_E,X_F)
        h_B = h_sw(T_E,P_E,X_B)
        h_D = CP.PropsSI('H','T',T_E,'Q',1,'Water')
        U_E = (1939.1 + 1.40562 * (T_S - 273.15) - 0.02075255 * (T_S - 273.15)**2 + 0.0023186 * (T_S - 273.15)**3)
        eq1 = F - (B + D)
        eq2 = F*X_F - B*X_B
        eq3 = m_s * hfg_S + F * h_F - D * h_D - B * h_B
        eq4 = m_s * hfg_S - U_E * A_E * (T_S - T_E)
        return [eq1,eq2,eq3,eq4]

    initial_guess_effect_2 = [D, B, X_B, T_E]

    sol_effect_2 = fsolve(Effect_mass_energyBalance,initial_guess_effect_2)

    D, B, X_B, T_E = sol_effect_2
    net_distillate += m_s
    T_E_arr.append(T_E)
    T_S_arr.append(T_S)
    F_arr.append(F)
    P_E_arr.append(CP.PropsSI('P','T',T_E,'Q',1,'Water')/1e6)
    m_s_arr.append(m_s)

    print(f"Effect {N_E}")
    print("************")
    print("m_s = ",m_s)
    print("D = ", D)
    print("B = ", B)
    print("X_B = ", X_B)
    print("T_E = ", T_E-273.15)
    print("BPE = ", BPE_sw(T_E, X_B), " K")
    print("Net distillate produced = ",net_distillate)
    print("")

print("Check final mass balance")
print("B_n + D_n + total distillate == F + m_s ---> ", (F1 + m_s_1 - B - D - net_distillate))
print("Productivity = ", (net_distillate - m_s_1)/F1)
# plt.plot([i for i in range(1,9)],T_E_arr,marker='o',c='k')
# plt.plot([i for i in range(1,9)],T_S_arr,marker='o',c='r')
plt.plot([i for i in range(1,9)],F_arr,marker='d',c='r')
plt.plot([i for i in range(1,9)],m_s_arr,marker='^',c='magenta')
# plt.plot([i for i in range(1,9)],P_E_arr,marker='x',c='b')
plt.show()
