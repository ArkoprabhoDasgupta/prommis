from pyomo.environ import ConcreteModel, SolverFactory, Var, Binary, Param
from idaes.core.solvers import get_solver
import numpy as np


def D_calculation(element_name, dosage ): #, pH_value):

    Element = ["Y", "Nd", "Dy", "Sm", "Gd", "Ce"]
    #System = ["5% dehpa 10% tbp", "2% dehpa 10% tbp", "2% dehpa", "2% cyanex"]

    m = ConcreteModel()

    m.a = Var(Element, domain=Binary)
    m.M = Var()
    m.B = Var()
    m.D = Var()

    slope = np.array([2.39, 0.5, 2.03, 0.7, 1.12, 0.55])

    intercept = np.array([-2.32, -1.89, -2.44, -2, -2.22, -1.94])

    LR_parameters = {}

    for i_e, element in enumerate(Element):
        LR_parameters[element] = (slope[i_e], intercept[i_e])

    @m.Constraint()
    def specify_element(m):
        return m.a[element_name] == 1

    @m.Constraint()
    def sum_of_binary_terms(m):
        return sum(m.a[e] for e in Element) == 1

    @m.Constraint()
    def slope_calculation(m):
        return m.M == sum(m.a[e] * LR_parameters[e][0] for e in Element)

    @m.Constraint()
    def intercept_calculation(m):
        return m.B == sum(m.a[e] * LR_parameters[e][1] for e in Element) + sum(
            m.a[e] * LR_parameters[e][0] for e in Element
        ) * np.log10((dosage + 10) / 15)

    # @m.Constraint()
    # def D_calculation(m):
    #     return m.D == 10 ** (m.M * pH_value + m.B)

    solver = SolverFactory("ipopt")
    solver.solve(m)

    return m.M(), m.B()
