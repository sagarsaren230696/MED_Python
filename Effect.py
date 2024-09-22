from pyomo.environ import *
from pyomo.core import *

import numpy as np

m = ConcreteModel()

m.Effect = Block()

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

m.P0 = Param(initialize=0.101) # reference pressure in MPa

## Mass flow rates
m.Effect.F = Var(initialize=2)
m.Effect.B = Var(initialize=1)
m.Effect.D = Var(initialize=1)

## Salinity
m.Effect.XF = Var(initialize=0.042)
m.Effect.XB = Var(initialize=0.042)

## Temperature
m.Effect.Ts = Var(initialize=353.15)
m.Effect.TF = Var(initialize=333.15)
m.Effect.TE = Var(initialize=343.15,bounds=(283.15,393.15))

## Pressure 
m.Effect.P = Var(initialize=0.101) # in MPa

## Enthalpy
m.Effect.hD = Var(initialize=293036)
m.Effect.hF_sw = Var(initialize=293036)
m.Effect.hF_sw0 = Var(initialize=293036)
m.Effect.hF_w = Var(initialize=293036)
m.Effect.hB_sw = Var(initialize=293036)
m.Effect.hB_sw0 = Var(initialize=293036)
m.Effect.hB_w = Var(initialize=293036)

## Property calculation
m.Effect.hD_eqn = Constraint(expr=m.Effect.hD == 141.355 + 4202.07 * (m.Effect.TE - 273.15) - 0.535 * (m.Effect.TE - 273.15)**2 + 0.004 * (m.Effect.TE - 273.15)**3)
## Feed water enthalpy
m.Effect.hF_w_eqn = Constraint(expr=m.Effect.hF_w == 141.355 + 4202.07 * (m.Effect.TF - 273.15) - 0.535 * (m.Effect.TF - 273.15)**2 + 0.004 * (m.Effect.TF - 273.15)**3)
m.Effect.hF_sw0_eqn = Constraint(expr=m.Effect.hF_sw0 == m.Effect.hF_w - m.Effect.XF * (
    -b[0] + b[1]*m.Effect.XF + b[2]*m.Effect.XF**2 + b[3]*m.Effect.XF**3 
    + b[4] * (m.Effect.TF - 273.15) + b[5] * (m.Effect.TF - 273.15)**2 + b[6] * (m.Effect.TF - 273.15) ** 3
    + b[7] * m.Effect.XF * (m.Effect.TF - 273.15) + b[8] * m.Effect.XF**2 * (m.Effect.TF - 273.15) + b[9] * m.Effect.XF * (m.Effect.TF - 273.15)**2 
))
m.Effect.hF_sw_eqn = Constraint(expr=m.Effect.hF_sw == m.Effect.hF_sw0 + (m.Effect.P - m.P0) * 
                                (a[0] + a[1] * (m.Effect.TF - 273.15) + a[2] * (m.Effect.TF - 273.15)**2 + a[3] * (m.Effect.TF - 273.15)**3 
                                 + m.Effect.XF * 1000 * (a[4] + a[5] * (m.Effect.TF - 273.15) + a[6] * (m.Effect.TF - 273.15)**2 + a[7] * (m.Effect.TF - 273.15)**3)))
## Brine water enthalpy
m.Effect.hB_w_eqn = Constraint(expr=m.Effect.hB_w == 141.355 + 4202.07 * (m.Effect.TE - 273.15) - 0.535 * (m.Effect.TE - 273.15)**2 + 0.004 * (m.Effect.TE - 273.15)**3)
m.Effect.hB_sw0_eqn = Constraint(expr=m.Effect.hB_sw0 == m.Effect.hB_w - m.Effect.XF * (
    -b[0] + b[1]*m.Effect.XF + b[2]*m.Effect.XF**2 + b[3]*m.Effect.XF**3 
    + b[4] * (m.Effect.TE - 273.15) + b[5] * (m.Effect.TE - 273.15)**2 + b[6] * (m.Effect.TE - 273.15) ** 3
    + b[7] * m.Effect.XF * (m.Effect.TE - 273.15) + b[8] * m.Effect.XF**2 * (m.Effect.TE - 273.15) + b[9] * m.Effect.XF * (m.Effect.TE - 273.15)**2 
))
m.Effect.hB_sw_eqn = Constraint(expr=m.Effect.hB_sw == m.Effect.hB_sw0 + (m.Effect.P - m.P0) * 
                                (a[0] + a[1] * (m.Effect.TE - 273.15) + a[2] * (m.Effect.TE - 273.15)**2 + a[3] * (m.Effect.TE - 273.15)**3 
                                 + m.Effect.XF * 1000 * (a[4] + a[5] * (m.Effect.TE - 273.15) + a[6] * (m.Effect.TE - 273.15)**2 + a[7] * (m.Effect.TE - 273.15)**3)))

## Mass and energy balance equations
m.Effect_massBalance = Constraint(expr=m.Effect.F == m.Effect.B + m.Effect.D)
m.Effect_salinityBalance = Constraint(expr=m.Effect.F*m.Effect.XF == m.Effect.B*m.Effect.XB)

## Fix the variable values
m.Effect.F.fix(2.0)
m.Effect.B.fix(1.0)
m.Effect.XF.fix(0.042)
m.Effect.TE.fix(343.15)
m.Effect.TF.fix(333.15)
m.Effect.P.fix(0.2)

solver = SolverFactory('ipopt')
solver.solve(m)

print("F = ",value(m.Effect.F))
print("B = ",value(m.Effect.B))
print("D = ",value(m.Effect.D))
print("XF = ",value(m.Effect.XF))
print("XB = ",value(m.Effect.XB))
print("hD = ",value(m.Effect.hD), " J/kg")
print("hF_w = ",value(m.Effect.hF_w), " J/kg")
print("hF_sw0 = ",value(m.Effect.hF_sw0), " J/kg")
print("hF_sw = ",value(m.Effect.hF_sw), " J/kg")