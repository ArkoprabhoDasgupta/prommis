"""
This is a Pyomo model which is used to describe the solvent extraction process of 
5 rare earth elements Y, Nd, Dy, Gd, Sm.

This solvent extraction process consists of 1 loading stage operated at a pH of 1.542 
and 1 stripping stage operated at pH 0.241. This operation is performed for the trial
experiment described in the phase 1 report, for a system of 5% DEHPA and 10% TBP extractant
dosage.

It is assumed that the volumes of the two phases are equal to 50 ml and there is no change 
in their amounts during the entire process.

In the extraction operation, the barren organic phase is put in contact with the feed aqueous
phase, and the rare earth elements are extracted into the organic phase. 
Then the loaded organic phase is used as a feed for the stripping operation, with a barren 
aqueous feed at the required pH, and the rare earth elements are stripped back in the aqueous
phase.

"""

from pyomo.environ import ConcreteModel, SolverFactory, Var
from idaes.core.solvers import get_solver

# Elements
Elements = ["Y", "Dy", "Gd", "Nd", "Sm"]

# pH correlation parameters
slope_pH = [2.39, 2.03, 1.12, 0.5, 0.7]
intercept_pH = [-2.32, -2.44, -2.22, -1.89, -2]
a = dict(zip(Elements, slope_pH))
b = dict(zip(Elements, intercept_pH))

# REE concentration in aqueous feed loading in ppm
C_loading_feed = [225, 60, 80, 75, 50]
C_aq_feed = dict(zip(Elements, C_loading_feed))

# pH values at loading and stripping
pH_loading = 1.5421
pH_stripping = 0.241

# Volumes of the phases in mL
V_aq = 50 
V_org = 50

m = ConcreteModel()
m.C_loading_aq = Var(Elements, bounds=(0,None))
m.C_loading_org = Var(Elements, bounds=(0,None))
m.C_stripping_aq = Var(Elements, bounds=(0,None))
m.C_stripping_org = Var(Elements, bounds=(0,None))
m.C_stripping_feed = Var(Elements, bounds=(0,None))
m.D_loading = Var(Elements, bounds=(0,None))
m.D_stripping = Var(Elements, bounds=(0,None))
m.E = Var(Elements, bounds=(0,100))
m.S = Var(Elements, bounds=(0,100))

# Loading

# D constraint loading
@m.Constraint(Elements)
def D_loading_ratio(m, e):
    return m.C_loading_org[e] == m.D_loading[e]*m.C_loading_aq[e]

# D calculation loading
@m.Constraint(Elements)
def D_loading_cal(m, e):
    return m.D_loading[e] == 10**(a[e]*pH_loading + b[e])

# Material balance loading
@m.Constraint(Elements)
def material_balance_load(m, e):
    return V_aq*C_aq_feed[e] == V_aq*m.C_loading_aq[e] + V_org*m.C_loading_org[e]

# Extraction percentage
@m.Constraint(Elements)
def extraction_load(m, e):
    return m.C_loading_aq[e]/C_aq_feed[e] == 1 - m.E[e]/100

# Stripping feed
@m.Constraint(Elements)
def stripping_feed(m, e):
    return m.C_stripping_feed[e]  == m.C_loading_org[e]

# D constraint stripping
@m.Constraint(Elements)
def D_stripping_ratio(m, e):
    return m.C_stripping_org[e] == m.D_stripping[e]*m.C_stripping_aq[e]

# D calculation stripping
@m.Constraint(Elements)
def D_stripping_cal(m, e):
    return m.D_stripping[e] == 10**(a[e]*pH_stripping + b[e])

# Material balance stripping
@m.Constraint(Elements)
def material_balance_strip(m, e):
    return V_org*m.C_stripping_feed[e] == V_aq*m.C_stripping_aq[e] + V_org*m.C_stripping_org[e]

# Stripping percentage
@m.Constraint(Elements)
def stripping_load(m, e):
    return m.C_stripping_org[e]/m.C_stripping_feed[e] == 1 - m.S[e]/100

solver = SolverFactory("ipopt")
#solver = get_solver("ipopt")
solver.solve(m, tee=True)

for e in Elements:
    print('Extraction percentage of ' + e + ' = ' + str(m.E[e]()) + ' %')

for e in Elements:
    print('Stripping percentage of ' + e + ' = ' + str(m.S[e]()) + ' %')