#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import (
    ConcreteModel,
    units,
    TransformationFactory,
    RangeSet,
    Var,
    minimize,
)
import pandas as pd
from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)
from idaes.core.util import to_json
import matplotlib.pyplot as plt
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution_new import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction

from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
    MixerSettlerExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package_for_optimizing import (
    SolventExtractionReactions,
)
from sklearn.metrics import r2_score

df = pd.read_excel("data for parmest.xlsx", sheet_name="Sheet2")


def build_model():
    """
    Method to build a steady state model for mixer settler solvent extraction
    Args:
        dosage: percentage dosage of extractant to the system.
        number_of_stages: number of stages in the mixer settler model.
    Returns:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    """

    m = ConcreteModel()

    m.system = RangeSet(0, 17)

    m.fs = FlowsheetBlock(m.system, dynamic=False)

    m.prop_o = REESolExOgParameters()
    m.leach_soln = LeachSolutionParameters()
    # m.leach_soln._has_inherent_reactions = False
    m.reaxn = SolventExtractionReactions()

    for s in m.system:

        if s in RangeSet(0, 14):

            m.fs[s].solex = SolventExtraction(
                number_of_finite_elements=1,
                aqueous_stream={
                    "property_package": m.leach_soln,
                    "flow_direction": FlowDirection.forward,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
                organic_stream={
                    "property_package": m.prop_o,
                    "flow_direction": FlowDirection.backward,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
                heterogeneous_reaction_package=m.reaxn,
                has_holdup=True,
            )

        else:

            m.fs[s].solex = SolventExtraction(
                number_of_finite_elements=2,
                aqueous_stream={
                    "property_package": m.leach_soln,
                    "flow_direction": FlowDirection.forward,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
                organic_stream={
                    "property_package": m.prop_o,
                    "flow_direction": FlowDirection.backward,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
                heterogeneous_reaction_package=m.reaxn,
                has_holdup=True,
            )

    return m


def set_inputs(m, df):
    """
    Set inlet conditions to the mixer settler solvent extraction model and fixing
    the parameters of the model.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
        dosage: percentage dosage of extractant to the system.
    Returns:
        None

    """

    for s in m.system:
        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(
            10 ** -df.loc[s, "pH"] * 1e3
        )
        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(
            10 ** -df.loc[s, "pH"] * 35e3
        )
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
            975.8e3 * df.loc[s, "dosage"] / 100
        )
        # m.fs.mixer_settler_ex[s].organic_inlet.extractant_dosage.fix(dosage_set[s])

        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(1e-7)
        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-7)

        m.fs[s].solex.aqueous_inlet.flow_vol.fix(62.01)

        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)

        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Al_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Ca_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Sc_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Y_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "La_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Ce_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Pr_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Nd_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Sm_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Gd_o"].fix(1e-7)
        m.fs[s].solex.organic_inlet.conc_mass_comp[0, "Dy_o"].fix(1e-7)

        m.fs[s].solex.organic_inlet.flow_vol.fix(62.01)

        m.fs[s].solex.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
        m.fs[s].solex.mscontactor.aqueous_inlet_state[:].temperature.fix(
            305.15 * units.K
        )
        m.fs[s].solex.mscontactor.organic[:, :].temperature.fix(305.15 * units.K)
        m.fs[s].solex.mscontactor.organic_inlet_state[:].temperature.fix(
            305.15 * units.K
        )

        m.fs[s].solex.mscontactor.volume[:].fix(0.4 * units.m**3)
        m.fs[s].solex.area_cross_stage[:] = 1
        m.fs[s].solex.elevation[:] = 0

        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Al"].fix(1e-7)
        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(1e-7)
        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(1e-7)
        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(1e-7)

        if s in RangeSet(0, 14):

            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(225)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(5)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(50)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(12.5)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(75)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(45)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(80)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(60)

        else:
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.339)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(2.082)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(4.982)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.737)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(2.093)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.257)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.561)
            m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.093)

    for e in ["Al", "Sc", "Ca", "Fe"]:
        m.reaxn.B0[e].fix(0)
        m.reaxn.B1[e].fix(0)
        m.reaxn.m0[e].fix(0)
        m.reaxn.m1[e].fix(0)

    m.fs[:].solex.distribution_extent_constraint[:, :, e].deactivate()
    m.fs[:].solex.mscontactor.heterogeneous_reaction_extent[
        0.0, :, f"{e}_mass_transfer"
    ].fix(0)


m = build_model()
set_inputs(m, df)

print(degrees_of_freedom(m))


REE_list = [
    e
    for e in m.leach_soln.component_list
    if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
]


@m.Expression(m.system, REE_list)
def percentage_extraction(m, s, e):
    return (
        (
            m.fs[s].solex.organic_outlet.conc_mass_comp[0, f"{e}_o"]
            * m.fs[s].solex.organic_outlet.flow_vol[0]
        )
        - (
            m.fs[s].solex.organic_inlet.conc_mass_comp[0, f"{e}_o"]
            * m.fs[s].solex.organic_inlet.flow_vol[0]
        )
    ) / (
        m.fs[s].solex.aqueous_inlet.conc_mass_comp[0, e]
        * m.fs[s].solex.aqueous_inlet.flow_vol[0]
    )


# @m.Objective(sense=minimize)
# def mse(m):

#     return sum(
#         (m.percentage_extraction[s, e] - df.loc[s, f"E {e}"] / 100) ** 2
#         for s in m.system
#         for e in REE_list
#         if df.loc[s, f"w_{e}"] == 1
#     ) / sum(
#         df.loc[s, f"w_{e}"]
#         for s in m.system
#         for e in REE_list
#         if df.loc[s, f"w_{e}"] == 1
#     )


@m.Objective(sense=minimize)
def mse(m):

    return sum(
        df.loc[s, f"w_{e}"]
        * (m.percentage_extraction[s, e] - df.loc[s, f"E {e}"] / 100) ** 2
        for s in m.system
        for e in REE_list
        if df.loc[s, f"w_{e}"] != 0
    ) / sum(df.loc[s, f"w_{e}"] for s in m.system for e in REE_list)


print(degrees_of_freedom(m))

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 5000
solver.options["halt_on_ampl_error"] = "yes"
solver.solve(m, tee=True)

for e in REE_list:
    Y_exp = [df.loc[s, f"E {e}"] for s in m.system if df.loc[s, f"w_{e}"] != 0]
    Y_model = [
        m.percentage_extraction[s, e]() * 100
        for s in m.system
        if df.loc[s, f"w_{e}"] != 0
    ]

    plt.plot(Y_exp, Y_exp)
    plt.scatter(Y_exp, Y_model)
    plt.title(f"{e}, {r2_score(Y_exp, Y_model)}")
    plt.xlabel("experiment")
    plt.ylabel("model")
    plt.show()
