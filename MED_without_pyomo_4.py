import numpy as np
from scipy.optimize import fsolve,minimize
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
    
def MED(args,N_effects,plotShow=False,printData=False):
    """
    args: F, D, A_E, T_S, T_F
    """
    F,D,A_E,T_S,T_F = args
    # print("Guesses : D = ",D," A_E = ",A_E)
    
    # F = 2.0
    # Feed flow rate in Effect 1
    F1 = F                                                  
    # D = 0.171
    # Feed salinity in Effect 1
    X_F = 0.042                                             
    X_F_1 = X_F                                             
    # T_E = 78 + 273.15
    # T_F = 333.15 #66.81 + 273.15
    # T_S = 353.15

    # Mass balance in Effect 1
    def Effect_massBalance_1(x):                            
        B, X_B = x
        eq1 = F - (B + D)
        eq2 = F*X_F - B*X_B
        return [eq1, eq2]

    # Initial guesses for B and X_B in Effect 1
    initial_guess = [1.0, 0.05]                             

    # Solving the mass balance equation in Effect 1
    sol = fsolve(Effect_massBalance_1,initial_guess)        
    
    # Getting B and X_B for Effect 1
    B, X_B = sol                                            

    # Specific enthalpy of vaporization of water at temperature T_S
    hfg_S = CP.PropsSI('H','T',T_S,'Q',1,'Water') - CP.PropsSI('H','T',T_S,'Q',0,'Water')

    # print("Effect pressure is ", P_E)

    # Total heat transfer coefficient at T_S
    U_E = (1939.1 + 1.40562 * (T_S - 273.15) - 0.02075255 * (T_S - 273.15)**2 + 0.0023186 * (T_S - 273.15)**3)
    # A_E = 90

    # Energy balance equation in Effect 1
    def Effect_energyBalance_1(x):
        m_s,T_E = x
        P_E = CP.PropsSI('P','T',T_E,'Q',1,'Water') / 1e6
        h_D = CP.PropsSI('H','T',T_E,'Q',1,'Water')
        h_F = h_sw(T_F,P_E,X_F)
        h_B = h_sw(T_E,P_E,X_B)
        eq1 = m_s * hfg_S + F * h_F - D * h_D - B * h_B
        eq2 = m_s*hfg_S - A_E * U_E * (T_S - T_E)
        return [eq1, eq2]

    # Initial guesses for m_s and T_E in Energy balance equation 
    initial_guess_energy = (2.0,78.0+273.15)

    # Solving energy balance equation 
    sol_energy = fsolve(Effect_energyBalance_1,initial_guess_energy)

    # Getting the m_s and T_E values
    m_s,T_E = sol_energy

    # print("Mass flow rate of steam is ", m_s, " kg/s" )
    # print("T_E = ",T_E-273.15)
    
    # just storing the value for Effect 1
    m_s_1 = m_s 

    # Pressure in Effect 1
    P_E = CP.PropsSI('P','T',T_E,'Q',1,'Water') / 1e6

    # Heat released from steam in Effect 1
    Q_E_1 = m_s_1 * hfg_S


    #######################################################
    ######### Effect 2 onwards ----------------------------
    #######################################################
    # N_effects = 8
    net_distillate = m_s_1 # net distillate amount till now 
    T_E_arr = [T_E]         # initializing Effect temperature array across Effects
    T_S_arr = [T_S]         # initializing Steam temperature array across Effects
    F_arr = [F]             # initializing Feed flow rate array across Effects
    P_E_arr = [P_E]         # initializing Effect pressure array across Effects
    m_s_arr = [m_s]         # initializing mass flow rate of steam array across Effects
    Q_E_arr = [Q_E_1]       # initializing heat transferred amount array across Effects

    if printData:
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

    # Looping over 2nd to Nth Effect
    for N_E in range(2,N_effects+1): 
        F = abs(B)                      # Brine outlet of previous effect is the feed to next effect
        m_s = D                         # Distillate flow rate is the steam mass flow rate in next effect
        T_S = T_E - BPE_sw(T_E, X_B)    # Decreasing in next Effect's steam temperature by the BPE temperature calculated at the previous effect temperature and brine salinity
        T_F = T_E                       # Feed temperature in next Effect is the temperature of the previous Effect
        X_F = X_B                       # Feed salinity of the next Effect is the brine salinity of the previous Effect

        # h_fg value at the temperature T_S
        hfg_S = CP.PropsSI('H','T',T_E,'Q',1,'Water') - CP.PropsSI('H','T',T_S,'Q',0,'Water')

        # Mass and energy balance equations
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

        # Initial guesses for the mass and energy balance equations
        initial_guess_effect_2 = [D, B, X_B, T_E]

        # Solving the mass and energy balance equations
        sol_effect_2 = fsolve(Effect_mass_energyBalance,initial_guess_effect_2)

        # Extracting the results for distillate flow rate, brine flow rate, brine salinity, and Effect temperature
        D, B, X_B, T_E = sol_effect_2
        net_distillate += m_s               # Adding the condensed m_s to the net distillate amount
        T_E_arr.append(T_E)                 # Appending the current Effect temperature to the T_E array
        T_S_arr.append(T_S)                 # Appending the current steam temperature to the T_S array
        F_arr.append(F)                     # Appending the feed water flow rate to the F array
        P_E_arr.append(CP.PropsSI('P','T',T_E,'Q',1,'Water')/1e6) # Appending the Effect pressure to the P_E array
        m_s_arr.append(m_s)                 # Appending the steam mass flow rate to the m_s array
        Q_E_arr.append(m_s * hfg_S)         # Appending the energy transferred value to the Q_E array 

        if printData:
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
    
    ### Add condenser
    T_sw_in = 25 + 273.15                   # Condenser seawater inlet temperature 
    T_sw_out = T_sw_in + 40 # K             # Condenser seawater outlet temperature
    T_S = T_E - BPE_sw(T_E, X_B)            # Steam temperature reduced by the BPE 
    hfg_S = CP.PropsSI('H','T',T_E,'Q',1,'Water') - CP.PropsSI('H','T',T_S,'Q',0,'Water')
    h_sw_in = h_sw(T_sw_in,0.101,X_F_1)     # Enthalpy value of seawater inlet
    h_sw_out = h_sw(T_sw_out,0.101,X_F_1)   # Enthalpy value of seawater outlet
    # print(D * hfg_S  , (h_sw_out - h_sw_in))
    F_sw = D * hfg_S / (h_sw_out - h_sw_in) # Required seawater feed flow rate for complete condensation of distillate 
    
    # print("Check final mass balance")
    # print("B_n + D_n + total distillate == F + m_s ---> ", (F1 + m_s_1 - B - D - net_distillate))
    # print("Productivity = ", (net_distillate - m_s_1)/F1)
    if plotShow:
        plt.plot([i for i in range(1,9)],T_E_arr,marker='o',c='k')
        plt.plot([i for i in range(1,9)],T_S_arr,marker='o',c='r')
        # plt.plot([i for i in range(1,9)],F_arr,marker='d',c='r')
        # plt.plot([i for i in range(1,9)],m_s_arr,marker='^',c='magenta')
        # plt.plot([i for i in range(1,9)],P_E_arr,marker='x',c='b')
        plt.show()

    Productivity = (net_distillate - m_s_1)/F1 # Productivity 

    PR = sum(Q_E_arr[1:])/Q_E_1                # Performance ratio
    
    return [Productivity,PR,F_sw]           

# Productivity,PR,F_sw = MED([2.0,0.171,90,353.15,333.15],8)
# print(Productivity,PR,F_sw)

# for _ in range(20):
#     Productivity,PR,F_sw = MED([F_sw,0.171,90,353.15,333.15],8)
#     print(Productivity,PR,F_sw)

# Productivity,PR,F_sw = MED([F_sw,0.171,90,353.15,333.15],8,plotShow=True)
# print(Productivity,PR,F_sw)

def obj(x,N_effects):
    prod,pr,f_sw = MED(x,N_effects=N_effects,plotShow=False)
    return -pr

x0 = [2.0,0.171,90,353.15,333.15]

prod_arr = []
for n in range(2,9):
    res = minimize(obj,x0,args=(n,),method="nelder-mead", options={'xatol': 1e-8, 'disp': True},bounds=((2.0,3.0),(0.12,0.24),(60,100),(343.15,363.15),(323.15,343.15)))

    print("nelder-mead : ",res.x)

    prod_arr.append(MED(res.x,n)[1])

print(prod_arr)
plt.scatter([n for n in range(2,9)],prod_arr)
plt.xlim([1,10])
plt.ylim([1,10])
plt.show()


# prod_arr = []

# T_arr = np.linspace(333.15,373.15,11)
# for T in T_arr:
#     Productivity,PR = MED([2.0,0.171,90,T,333.15],8,plotShow=False)

#     prod_arr.append(Productivity * 2.0)

# # print(prod_arr)
# plt.scatter(T_arr-273.15,prod_arr)
# plt.show()
# res = minimize(obj,x0,method="SLSQP", options={'xatol': 1e-8, 'disp': True},bounds=((2.0,3.0),(0.12,0.24),(60,100),(343.15,363.15)))

# print("SLSQP : ",res.x)

# res = minimize(obj,x0,method="L-BFGS-B", options={'xatol': 1e-8, 'disp': True},bounds=((2.0,3.0),(0.12,0.24),(60,100),(343.15,363.15)))

# print("L-BFGS-B : ",res.x)
# print(MED(res.x,plotShow=True,printData=False))

