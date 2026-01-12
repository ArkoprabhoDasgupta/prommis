#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import ConcreteModel, units, Set, TransformationFactory
import matplotlib.pyplot as plt
from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)
from idaes.core.util import to_json

from idaes.core.util.scaling import set_scaling_factor
from idaes.core.solvers import get_solver

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution_new import REESolExOgParameters
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
    MixerSettlerExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package_new_modified import (
    SolventExtractionReactions,
)


"""
Method to build a steady state model for mixer settler solvent extraction
Args:
    dosage: percentage dosage of extractant to the system.
    number_of_stages: number of stages in the mixer settler model.
Returns:
    m: ConcreteModel object with the mixer-settler solvent extraction system.
"""

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()

# m.fs.reaxn.extractant_dosage = dosage

number_of_stages = 1

m.scenario = Set(initialize=["P1D1", "P2D1"])

m.fs.mixer_settler_ex = MixerSettlerExtraction(
    m.scenario,
    number_of_stages=number_of_stages,
    aqueous_stream={
        "property_package": m.fs.leach_soln,
        "flow_direction": FlowDirection.forward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
    organic_stream={
        "property_package": m.fs.prop_o,
        "flow_direction": FlowDirection.backward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
    heterogeneous_reaction_package=m.fs.reaxn,
    has_holdup=True,
    settler_transformation_method="dae.finite_difference",
    settler_transformation_scheme="BACKWARD",
    settler_finite_elements=4,
)


"""
Set inlet conditions to the mixer settler solvent extraction model and fixing
the parameters of the model.
Args:
    m: ConcreteModel object with the mixer-settler solvent extraction system.
    dosage: percentage dosage of extractant to the system.
Returns:
    None

"""
pH_set = {
    "P1D1": 1.3,
    "P2D1": 1.9,
}

for s in m.scenario:
    m.fs.mixer_settler_ex[s].aqueous_inlet.conc_mass_comp[0, "H"].fix(
        10 ** -pH_set[s] * 1e3
    )
    m.fs.mixer_settler_ex[s].aqueous_inlet.conc_mass_comp[0, "SO4"].fix(
        10 ** -pH_set[s] * 48e3
    )

dosage = 25

m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-8)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Al"].fix(1e-8)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Ca"].fix(1e-8)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Fe"].fix(1e-7)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.277)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.346)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "La"].fix(2.091)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Ce"].fix(5.013)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.734)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Nd"].fix(2.106)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.235)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.567)
m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.096)

m.fs.mixer_settler_ex[:].aqueous_inlet.flow_vol.fix(60.01)

m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
    975.8e3 * dosage / 100
)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Al_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Ca_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Fe_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Sc_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Y_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "La_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Ce_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Pr_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Nd_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Sm_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Gd_o"].fix(1e-10)
m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Dy_o"].fix(1e-10)

m.fs.mixer_settler_ex[:].organic_inlet.flow_vol.fix(60.01)

m.fs.mixer_settler_ex[:].mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.mixer_settler_ex[:].mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
    305.15 * units.K
)

m.fs.mixer_settler_ex[:].mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)

m.fs.mixer_settler_ex[:].organic_settler[:].unit.area.fix(1)
m.fs.mixer_settler_ex[:].aqueous_settler[:].unit.area.fix(1)
m.fs.mixer_settler_ex[:].aqueous_settler[:].unit.length.fix(1)
m.fs.mixer_settler_ex[:].organic_settler[:].unit.length.fix(1)

# for e in ["Al", "Ca", "Fe", "Sc"]:
#     m.fs.mixer_settler_ex[:].mixer[1].unit.mscontactor.heterogeneous_reaction_extent[
#         0.0, 1, f"{e}_mass_transfer"
#     ].fix(0)
#     m.fs.mixer_settler_ex[:].mixer[1].unit.distribution_extent_constraint[
#         0, 1, e
#     ].deactivate()

# for s in m.scenario:
#     MixerSettlerExtractionInitializer().initialize(m.fs.mixer_settler_ex[s])

# for s in m.scenario:
#     for e in ["La", "Ce", "Nd"]:
#         # set_scaling_factor(
#         #     m.fs.mixer_settler_ex[8]
#         #     .mixer[1]
#         #     .unit.mscontactor.heterogeneous_reactions[0.0, 1]
#         #     .distribution_expression_constraint[e],
#         #     1e-1,
#         # )

#         set_scaling_factor(
#             m.fs.mixer_settler_ex[s]
#             .mixer[1]
#             .unit.distribution_extent_constraint[0.0, 1, e],
#             1e-2,
#         )

#         # set_scaling_factor(
#         #     m.fs.mixer_settler_ex[s]
#         #     .mixer[2]
#         #     .unit.distribution_extent_constraint[0.0, 1, e],
#         #     1e1,
#         # )

#         set_scaling_factor(
#             m.fs.mixer_settler_ex[s]
#             .mixer[1]
#             .unit.mscontactor.aqueous[0.0, 1]
#             .hso4_dissociation,
#             1e-2,
#         )

#         set_scaling_factor(
#             m.fs.mixer_settler_ex[s]
#             .mixer[2]
#             .unit.mscontactor.aqueous[0.0, 1]
#             .hso4_dissociation,
#             1e-2,
#         )

#         set_scaling_factor(
#             m.fs.mixer_settler_ex[s]
#             .mixer[1]
#             .unit.mscontactor.aqueous[0.0, 1]
#             .pH_constraint["liquid"],
#             1e2,
#         )

#         set_scaling_factor(
#             m.fs.mixer_settler_ex[s]
#             .mixer[2]
#             .unit.mscontactor.aqueous[0.0, 1]
#             .pH_constraint["liquid"],
#             1e2,
#         )

# scaling = TransformationFactory("core.scale_model")
# scaled_model = scaling.create_using(m, rename=False)

for s in m.scenario:
    MixerSettlerExtractionInitializer().initialize(m.fs.mixer_settler_ex[s])

solver = get_solver("ipopt_v2")
results = solver.solve(m, tee=True)

# solve_model(scaled_model)

# # for s in m.system:
# #     MixerSettlerExtractionInitializer().initialize(m.fs.mixer_settler_ex[s])

# scaling.propagate_solution(scaled_model, m)

# solver = get_solver("ipopt_v2")
# results = solver.solve(m, tee=True)

# if __name__ == "__main__":
#     m, results = main(dosage, number_of_stages)

percentage_extraction = {}

# for s in m.scenario:
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]:
        percentage_extraction[e] = [
            (
                (
                    m.fs.mixer_settler_ex[s].organic_outlet.conc_mass_comp[
                        0, f"{e}_o"
                    ]()
                    * m.fs.mixer_settler_ex[s].organic_outlet.flow_vol[0]()
                    - m.fs.mixer_settler_ex[s].organic_inlet.conc_mass_comp[
                        0, f"{e}_o"
                    ]()
                    * m.fs.mixer_settler_ex[s].organic_inlet.flow_vol[0]()
                )
                / (
                    m.fs.mixer_settler_ex[s].aqueous_inlet.conc_mass_comp[0, e]()
                    * m.fs.mixer_settler_ex[s].aqueous_inlet.flow_vol[0]()
                )
            )
            * 100
            for s in m.scenario
        ]

percentage_extraction["TREE"] = [
    (
        (
            sum(
                m.fs.mixer_settler_ex[s].organic_outlet.conc_mass_comp[0, f"{e}_o"]()
                for e in m.fs.leach_soln.component_list
                if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
            )
            * m.fs.mixer_settler_ex[s].organic_outlet.flow_vol[0]()
            - sum(
                m.fs.mixer_settler_ex[s].organic_inlet.conc_mass_comp[0, f"{e}_o"]()
                for e in m.fs.leach_soln.component_list
                if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
            )
            * m.fs.mixer_settler_ex[s].organic_inlet.flow_vol[0]()
        )
        / (
            sum(
                m.fs.mixer_settler_ex[s].aqueous_inlet.conc_mass_comp[0, e]()
                for e in m.fs.leach_soln.component_list
                if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
            )
            * m.fs.mixer_settler_ex[s].aqueous_inlet.flow_vol[0]()
        )
    )
    * 100
    for s in m.scenario
]

pH_list = [
    m.fs.mixer_settler_ex[s]
    .mixer[1]
    .unit.mscontactor.aqueous[0.0, 1]
    .pH_phase["liquid"]()
    for s in m.scenario
]

plt.scatter([1.5, 2, 2.5], [20, 35, 40], marker="o", color="red")
plt.plot(pH_list, percentage_extraction["La"], marker="v", color="black")
plt.xlabel("pH")
plt.ylabel("Extraction %")
plt.title(f"Extraction % comparison for 5% DEHPA 10% TBP")
plt.legend(fontsize=8)
