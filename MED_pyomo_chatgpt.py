# med_pyomo.py
# MED energy-balance simulation in Pyomo (algebraic, steady-state)

from pyomo.environ import (
    ConcreteModel, Var, Param, RangeSet, Constraint, NonNegativeReals, Reals,
    value, SolverFactory, Objective, minimize
)
from math import log

# -------------------------------
# Seawater enthalpy correlations (from your code, units as-is)
# -------------------------------
a = [
    996.7767, -3.2406, 0.0127, -4.7723e-5,
    -1.1748, 0.01169, -2.6185e-5, 7.0661e-8
]
b = [
    -2.34825e4, 3.15183e5, 2.80269e6, -1.44606e7,
    7.82607e3, -4.41733e1, 2.1394e-1, -1.99108e4,
    2.77846e4, 9.72801e1
]
c = [
    -4.5838530457E-04, 2.8230948284E-01, 1.7945189194E+01,
    1.5361752708E-04, 5.2669058133E-02, 6.5604855793E+00
]

def h_w(T):
    # pure water liquid enthalpy (J/kg), baseline 0 at 0°C-ish (your polynomial)
    return 141.355 + 4202.07 * (T - 273.15) - 0.535 * (T - 273.15)**2 + 0.004 * (T - 273.15)**3

def h_sw0(T, X):
    hw = h_w(T)
    DT = T - 273.15
    return hw - X * (
        -b[0] + b[1]*X + b[2]*X**2 + b[3]*X**3
        + b[4]*DT + b[5]*DT**2 + b[6]*DT**3
        + b[7]*X*DT + b[8]*X**2 * DT + b[9]*X*DT**2
    )

def h_sw(T, P, X):
    # P in MPa, T in K, X mass fraction salt
    DT = T - 273.15
    return h_sw0(T, X) + (P - 0.101) * (
        a[0] + a[1]*DT + a[2]*DT**2 + a[3]*DT**3
        + X*1000*(a[4] + a[5]*DT + a[6]*DT**2 + a[7]*DT**3)
    )

def BPE_sw(T, X):
    # Boiling point elevation (K)
    T_degC = T - 273.15
    A = c[0]*T_degC**2 + c[1]*T_degC + c[2]
    B = c[3]*T_degC**2 + c[4]*T_degC + c[5]
    return A*X**2 + B*X

# -------------------------------
# Steam/Water property surrogates (solver-friendly)
# -------------------------------
# Antoine for water (valid ~1–100°C); we’ll use it across range typical for MED.
# P_sat in kPa -> convert to MPa. T in K.
ANTOINE_A = 8.07131
ANTOINE_B = 1730.63
ANTOINE_C = 233.426

def Psat_MPaa_from_T(T):
    # returns MPa; clamp domain in caller if needed
    T_C = T - 273.15
    # Antoine uses mmHg; but these coefficients are for kPa-ish? To avoid mismatch,
    # we switch to a robust log form for kPa using Buck-like form:
    # Using simplified form: log10(P_mmHg) = A - B/(C + T_C)
    # Convert mmHg to kPa: 1 mmHg = 0.133322 kPa
    P_mmHg = 10.0**(ANTOINE_A - ANTOINE_B / (ANTOINE_C + max(T_C, 1.0)))
    P_kPa = P_mmHg * 0.133322
    return P_kPa / 1000.0  # MPa

# Watson correlation for latent heat:
# h_fg(T) = h_fg(T1) * ((1 - Tr)/(1 - Tr1))^0.38 ; Tc = 647.096 K, choose T1 = 373.15 K
T_c = 647.096
T1 = 373.15
h_fg_T1 = 2.257e6  # J/kg at 100°C

def h_fg(T):
    Tr = T / T_c
    Tr1 = T1 / T_c
    # Guard near Tc
    ratio = max(1e-6, (1.0 - Tr)) / (1.0 - Tr1)
    return h_fg_T1 * (ratio ** 0.38)

# Saturated vapor enthalpy ~ saturated liquid enthalpy + h_fg(T)
def h_g_sat(T):
    return h_w(T) + h_fg(T)

# Overall U correlation from your code (W/m2-K) as function of steam-side T_S
def U_E_of(T_S):
    DTc = T_S - 273.15
    return 1939.1 + 1.40562*DTc - 0.02075255*DTc**2 + 0.0023186*DTc**3

# -------------------------------
# Build and solve MED model
# -------------------------------
def build_med_model(
    N=8,
    F1=2.0,                 # feed (kg/s) to effect 1
    X_F1=0.042,             # mass fraction salt in feed
    T_S1=353.15,            # steam temperature to effect 1 (K)
    T_F1=333.15,            # feed brine temperature to effect 1 (K)
    A_E=90.0,               # area per effect (m2)
    T_sw_in=25+273.15,      # condenser seawater in (K)
    T_sw_out=None           # condenser seawater out (K), default T_in + 40 K
):
    if T_sw_out is None:
        T_sw_out = T_sw_in + 40.0

    m = ConcreteModel()
    m.I = RangeSet(1, N)

    # Parameters
    m.F1 = Param(initialize=F1)
    m.XF1 = Param(initialize=X_F1)
    m.TS1 = Param(initialize=T_S1)
    m.TF1 = Param(initialize=T_F1)
    m.AE  = Param(initialize=A_E)
    m.Tsw_in = Param(initialize=T_sw_in)
    m.Tsw_out = Param(initialize=T_sw_out)

    # Variables per effect
    m.D = Var(m.I, domain=NonNegativeReals, initialize=0.15)     # distillate (kg/s)
    m.B = Var(m.I, domain=NonNegativeReals, initialize=1.0)      # brine (kg/s)
    m.XB = Var(m.I, bounds=(1e-6, 0.25), initialize=0.05)        # brine salinity (mass fraction)
    m.TE = Var(m.I, bounds=(300.0, 390.0), initialize=351.0)     # effect boiling temperature (K)

    # Helper "link" variables (computed-like, but as Vars so equations remain symbolic)
    m.F  = Var(m.I, domain=NonNegativeReals)     # feed to effect i (kg/s)
    m.XF = Var(m.I, bounds=(1e-6, 0.25))         # feed salinity to effect i
    m.TS = Var(m.I, bounds=(300.0, 390.0))       # steam temp to effect i (K)
    m.ms = Var(m.I, domain=NonNegativeReals)     # condensing steam mass to effect i (kg/s)
    m.PE = Var(m.I, bounds=(0.001, 0.3))         # effect pressure (MPa)

    # First-effect feed/steam
    def first_links(m):
        yield m.F[1] == m.F1
        yield m.XF[1] == m.XF1
        yield m.TS[1] == m.TS1
    m.first_link1 = Constraint(expr=next(first_links(m)))
    m.first_link2 = Constraint(expr=next(first_links(m)))
    m.first_link3 = Constraint(expr=next(first_links(m)))

    # Subsequent-effect linking & steam source
    m.link_F = Constraint(m.I - {1}, rule=lambda m, i: m.F[i] == m.B[i-1])
    m.link_XF = Constraint(m.I - {1}, rule=lambda m, i: m.XF[i] == m.XB[i-1])
    m.link_ms = Constraint(m.I - {1}, rule=lambda m, i: m.ms[i] == m.D[i-1])
    # Steam temp is previous effect boiling minus BPE at previous brine
    def steam_temp_rule(m, i):
        if i == 1:
            return Constraint.Skip
        # TS[i] = TE[i-1] - BPE(TE[i-1], XB[i-1])
        return m.TS[i] == m.TE[i-1] - BPE_sw(m.TE[i-1], m.XB[i-1])
    m.link_TS = Constraint(m.I, rule=steam_temp_rule)

    # Pressure in each effect = Psat(TE)
    m.Psat = Constraint(m.I, rule=lambda m, i: m.PE[i] == Psat_MPaa_from_T(m.TE[i]))

    # Mass balances
    m.mb_total = Constraint(m.I, rule=lambda m, i: m.F[i] == m.B[i] + m.D[i])
    # Salt balance (distillate salt-free)
    m.mb_salt  = Constraint(m.I, rule=lambda m, i: m.F[i]*m.XF[i] == m.B[i]*m.XB[i])

    # Heat/energy balances
    # q_i = ms_i * h_fg(TS_i) = U*A*(TS_i - TE_i)
    m.ht = Constraint(m.I, rule=lambda m, i: m.ms[i]*h_fg(m.TS[i]) == U_E_of(m.TS[i]) * m.AE * (m.TS[i] - m.TE[i]))

    # Effect-level energy balance:
    # ms_i*h_fg(TS_i) + F_i*h_sw(TF_i, PE_i, XF_i) = D_i*h_g_sat(TE_i) + B_i*h_sw(TE_i, PE_i, XB_i)
    def TF_rule(m, i):
        # Feed temperature to effect i: TF1 for i=1, else previous TE
        return m.TF1 if i == 1 else m.TE[i-1]

    m.eb = Constraint(
        m.I,
        rule=lambda m, i: m.ms[i]*h_fg(m.TS[i]) + m.F[i]*h_sw(TF_rule(m, i), m.PE[i], m.XF[i]) ==
                          m.D[i]*h_g_sat(m.TE[i]) + m.B[i]*h_sw(m.TE[i], m.PE[i], m.XB[i])
    )

    # Initialize ms[1] loosely via UA*DT/hfg
    # (Not a constraint, just a better starting point.)
    m.ms[1].set_value(max(1e-3, U_E_of(value(m.TS1))*value(m.AE)*(value(m.TS1) - (value(m.TS1)-5.0)) / h_fg(value(m.TS1))))

    # Performance metrics
    # Q1 = ms1*hfg1 ; PR = sum_{i=2..N} (ms_i*hfg_i) / Q1 ; Productivity = (sum ms_i - ms1)/F1
    def Q_i(i):
        return m.ms[i]*h_fg(m.TS[i])

    m.Q1 = Var(initialize=1.0)
    m.Q1_def = Constraint(expr=m.Q1 == Q_i(1))

    # For convenience, compute PR & Productivity post-solve.

    # Optional: a very light objective to help solver (minimize TE_N just to have an objective)
    m.obj = Objective(expr=m.TE[N], sense=minimize)

    return m

def solve_med(
    N=8, F1=2.0, X_F1=0.042, T_S1=353.15, T_F1=333.15, A_E=90.0,
    T_sw_in=25+273.15, T_sw_out=None, solver='ipopt', tee=False
):
    m = build_med_model(N, F1, X_F1, T_S1, T_F1, A_E, T_sw_in, T_sw_out)
    opt = SolverFactory(solver)
    results = opt.solve(m, tee=tee)

    # Collect arrays
    D = [value(m.D[i]) for i in m.I]
    B = [value(m.B[i]) for i in m.I]
    XB = [value(m.XB[i]) for i in m.I]
    TE = [value(m.TE[i]) for i in m.I]
    TS = [value(m.TS[i]) for i in m.I]
    ms = [value(m.ms[i]) for i in m.I]
    PE = [value(m.PE[i]) for i in m.I]
    Qi = [ms[i-1]*h_fg(TS[i-1]) for i in m.I]

    net_distillate = sum(ms)  # matches your running sum of "m_s"
    Productivity = (net_distillate - ms[0]) / F1 if F1 > 0 else 0.0
    Q1 = Qi[0]
    PR = (sum(Qi[1:]) / Q1) if Q1 > 0 else 0.0

    # Simple condenser duty & cooling-water flow like your code (use last effect)
    T_E_last = TE[-1]
    X_B_last = XB[-1]
    T_S_last = T_E_last - BPE_sw(T_E_last, X_B_last)
    hfg_last = h_fg(T_S_last)
    h_sw_in = h_sw(T_sw_in, 0.101, X_F1)
    h_sw_out = h_sw(T_sw_out if T_sw_out is not None else T_sw_in + 40.0, 0.101, X_F1)
    F_sw = D[-1]*hfg_last / max(1e-6, (h_sw_out - h_sw_in))

    return {
        'status': str(results.solver.termination_condition),
        'D': D, 'B': B, 'XB': XB, 'TE': TE, 'TS': TS, 'ms': ms, 'PE': PE, 'Qi': Qi,
        'Productivity': Productivity,
        'PR': PR,
        'F_sw': F_sw
    }

if __name__ == "__main__":
    out = solve_med(
        N=8,
        F1=2.0,
        X_F1=0.042,
        T_S1=353.15,
        T_F1=333.15,
        A_E=90.0,
        T_sw_in=25+273.15
    )
    print("Status:", out['status'])
    print("PR =", out['PR'])
    print("Productivity =", out['Productivity'])
    print("Cooling seawater flow F_sw (kg/s) =", out['F_sw'])
    print("TE (°C):", [t-273.15 for t in out['TE']])
