from pyomo.environ import ConcreteModel, SolverFactory, Var, Binary, Param
from idaes.core.solvers import get_solver
import numpy as np


def D_calculation(element_name, system_name, pH_value):  # Add m to the arguments

    Element = ["Y", "Nd", "Dy", "Sm", "Gd", "Ce"]
    System = ["5% dehpa 10% tbp", "2% dehpa 10% tbp", "2% dehpa", "2% cyanex"]

    m = ConcreteModel()  # Remove m and do something like solex = m.fs.solex
    # Use solex.a, M, B etc
    m.a = Var(Element, System, domain=Binary)
    m.M = Var()
    m.B = Var()
    m.D = Var()

    slope = np.array(
        [
            [2.39, 0.5, 2.03, 0.7, 1.12, 0.55],
            [2.39, 0.5, 2.03, 0.7, 1.12, 0.55],
            [2.964, 0.894, 0.937, 1.711, 1.962, 0.512],
            [1.295, 0.269, 1.153, 0.492, 0.592, 0.278],
        ]
    )

    intercept = np.array(
        [
            [-2.32, -1.89, -2.44, -2, -2.22, -1.94],
            [-2.55, -1.93, -2.63, -2.06, -2.32, -1.99],
            [-1.16, -1.91, -1.45, -2.26, -2.12, -1.51],
            [-1.97, -1.44, -2.17, -1.85, -1.88, -1.48],
        ]
    )

    LR_parameters = {}

    for i_e, element in enumerate(Element):
        for i_s, system in enumerate(System):
            LR_parameters[element, system] = (slope[i_s, i_e], intercept[i_s, i_e])

    @m.Constraint()
    def specify_element(m):
        return m.a[element_name, system_name] == 1

    @m.Constraint()
    def sum_of_binary_terms(m):
        return sum(m.a[e, s] for e in Element for s in System) == 1

    @m.Constraint()
    def slope_calculation(m):
        return m.M == sum(
            m.a[e, s] * LR_parameters[e, s][0] for e in Element for s in System
        )

    @m.Constraint()
    def intercept_calculation(m):
        return m.B == sum(
            m.a[e, s] * LR_parameters[e, s][1] for e in Element for s in System
        )

    @m.Constraint()
    def D_calculation(m):
        return m.D == 10 ** (m.M * pH_value + m.B)

    # Remove solver stuff
    solver = SolverFactory("ipopt")
    solver.solve(m)

    return m.D()  # Return m only
