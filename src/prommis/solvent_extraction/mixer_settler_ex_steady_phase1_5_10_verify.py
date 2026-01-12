#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import ConcreteModel, units, TransformationFactory, RangeSet
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
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
    MixerSettlerExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package_new_modified import (
    SolventExtractionReactions,
)


def build_model(dosage, number_of_stages):
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
    m.fs.leach_soln._has_inherent_reactions = False
    m.fs.reaxn = SolventExtractionReactions()

    m.fs.reaxn.extractant_dosage = dosage

    m.system = RangeSet(1, 13)

    m.fs.mixer_settler_ex = MixerSettlerExtraction(
        m.system,
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

    return m


def set_inputs(m, dosage):
    """
    Set inlet conditions to the mixer settler solvent extraction model and fixing
    the parameters of the model.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
        dosage: percentage dosage of extractant to the system.
    Returns:
        None

    """
    pH_set = dict(
        zip(
            RangeSet(1, 13),
            sorted(
                [
                    0.25,
                    0.53,
                    0.81,
                    1.13,
                    1.4,
                    1.6,
                    1.68,
                    0.34,
                    0.7,
                    0.97,
                    1.26,
                    1.5,
                    1.64,
                ]
            ),
        )
    )

    # pH_set = {
    #     1: 0.25,
    #     2: 0.53,
    #     3: 0.81,
    #     4: 1.13,
    #     5: 1.4,
    #     6: 1.6,
    #     7: 1.68,
    # }
    dosage_set = {
        1: 5,
        2: 5,
        3: 5,
        4: 5,
        5: 5,
        6: 5,
        7: 5,
        8: 5,
        9: 5,
        10: 5,
        11: 5,
        12: 5,
        13: 5,
    }

    for s in m.system:
        m.fs.mixer_settler_ex[s].aqueous_inlet.conc_mass_comp[0, "H"].fix(
            10 ** -pH_set[s] * 1e3
        )
        m.fs.mixer_settler_ex[s].aqueous_inlet.conc_mass_comp[0, "Cl"].fix(
            10 ** -pH_set[s] * 35e3
        )
        m.fs.mixer_settler_ex[s].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
            975.8e3 * dosage_set[s] / 100
        )
        # m.fs.mixer_settler_ex[s].organic_inlet.extractant_dosage.fix(dosage_set[s])

    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
    # m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "H"].fix(10**-pH * 1e3)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "SO4"].fix(1e-7)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-7)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Al"].fix(1e-7)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Ca"].fix(1e-7)
    # m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Cl"].fix(10**-pH * 35e3)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Fe"].fix(1e-7)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Sc"].fix(1e-7)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Y"].fix(225)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "La"].fix(5)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Ce"].fix(50)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Pr"].fix(12.5)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Nd"].fix(75)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Sm"].fix(45)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Gd"].fix(80)
    m.fs.mixer_settler_ex[:].aqueous_inlet.conc_mass_comp[0, "Dy"].fix(60)

    m.fs.mixer_settler_ex[:].aqueous_inlet.flow_vol.fix(62.01)

    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)
    # m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
    #     975.8e3 * dosage / 100
    # )
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Al_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Ca_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Fe_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Sc_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Y_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "La_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Ce_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Pr_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Nd_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Sm_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Gd_o"].fix(1e-7)
    m.fs.mixer_settler_ex[:].organic_inlet.conc_mass_comp[0, "Dy_o"].fix(1e-7)

    m.fs.mixer_settler_ex[:].organic_inlet.flow_vol.fix(62.01)

    m.fs.mixer_settler_ex[:].mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
        305.15 * units.K
    )
    m.fs.mixer_settler_ex[:].mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
        305.15 * units.K
    )

    m.fs.mixer_settler_ex[:].mixer[:].unit.mscontactor.volume[:].fix(1e-4 * units.m**3)

    m.fs.mixer_settler_ex[:].organic_settler[:].unit.area.fix(1e-2)
    m.fs.mixer_settler_ex[:].aqueous_settler[:].unit.area.fix(1e-2)
    m.fs.mixer_settler_ex[:].aqueous_settler[:].unit.length.fix(1e-2)
    m.fs.mixer_settler_ex[:].organic_settler[:].unit.length.fix(1e-2)


def model_buildup_and_set_inputs(dosage, number_of_stages):
    """
    A function to build up the mixer settler solvent extraction model and set inlet
    streams to the model.
    Args:
        dosage: percentage dosage of extractant to the system.
        number_of_stages: Number of stages in the mixer settler model.
    Returns:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    """
    m = build_model(dosage, number_of_stages)
    set_inputs(m, dosage)

    return m


def initialize_steady_model(m):
    """
    A function to initialize the mixer settler solvent extraction model with the
    default initializer after setting input conditions.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None
    """
    initializer = m.fs.mixer_settler_ex.default_initializer()
    initializer.initialize(m.fs.mixer_settler_ex)


def solve_model(m):
    """
    A function to solve the initialized mixer-settler solvent extraction model.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None
    """
    solver = get_solver("ipopt_v2")
    results = solver.solve(m, tee=True)

    return results


def export_results(m):
    """
    A function to export the results of the solved mixer-settler solvent extraction
    model to a json file.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None
    """
    to_json(m, fname="mixer_settler_extraction.json", human_read=True)


def main(dosage, number_of_stages):
    """
    The main function used to build a mixer settler solvent extraction model, set inlets
    to the model, initialize the model, solve the model and export the results to a json
    file.
    Args:
        dosage: percentage dosage of extractant to the system.
        number_of_stages: Number of stages in the mixer settler model.
    Returns:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    """
    m = model_buildup_and_set_inputs(dosage, number_of_stages)
    initialize_steady_model(m)
    results = solve_model(m)
    export_results(m)

    return m, results


dosage = 5
number_of_stages = 1

m = model_buildup_and_set_inputs(dosage, number_of_stages)

m.fs.mixer_settler_ex[:].mixer[1].unit.mscontactor.aqueous[
    0.0, 1
].hso4_dissociation.deactivate()

for s in m.system:
    for z in m.fs.mixer_settler_ex[s].aqueous_settler[1].unit.length_domain:
        if z != 0:
            m.fs.mixer_settler_ex[s].aqueous_settler[1].unit.properties[
                0.0, z
            ].hso4_dissociation.deactivate()


for e in ["Al", "Ca", "Fe", "Sc"]:
    m.fs.mixer_settler_ex[:].mixer[1].unit.mscontactor.heterogeneous_reaction_extent[
        0.0, 1, f"{e}_mass_transfer"
    ].fix(0)
    m.fs.mixer_settler_ex[:].mixer[1].unit.distribution_extent_constraint[
        0, 1, e
    ].deactivate()


print(degrees_of_freedom(m))

# for s in m.system:
#     # set_scaling_factor(
#     #     m.fs.mixer_settler_ex[s]
#     #     .mixer[1]
#     #     .unit.mscontactor.aqueous_inlet_state[0.0]
#     #     .pH_constraint["liquid"],
#     #     1e-1,
#     # )

for s in m.system:
    for e in ["La", "Ce", "Pr", "Nd", "Sm","Gd","Dy","Y"]:
        # set_scaling_factor(
        #     m.fs.mixer_settler_ex[8]
        #     .mixer[1]
        #     .unit.mscontactor.heterogeneous_reactions[0.0, 1]
        #     .distribution_expression_constraint[e],
        #     1e-1,
        # )

        set_scaling_factor(
            m.fs.mixer_settler_ex[s]
            .mixer[1]
            .unit.distribution_extent_constraint[0.0, 1, e],
            1e1,
        )

# # set_scaling_factor(m.fs.mixer_settler_ex[:].mixer[1].unit.mscontactor.aqueous[0.0,1].conc_mol_comp['HSO4'])

scaling = TransformationFactory("core.scale_model")
scaled_model = scaling.create_using(m, rename=False)

# MixerSettlerExtractionInitializer().initialize(scaled_model.fs.mixer_settler_ex[8])

# for s in m.system:
#     MixerSettlerExtractionInitializer().initialize(scaled_model.fs.mixer_settler_ex[s])

solve_model(scaled_model)

# # initialize_steady_model(m)
# solve_model(scaled_model)

# for s in m.system:
#     MixerSettlerExtractionInitializer().initialize(m.fs.mixer_settler_ex[s])


scaling.propagate_solution(scaled_model, m)

# # if __name__ == "__main__":
# #     m, results = main(dosage, number_of_stages)

percentage_extraction = {}

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
            for s in m.system
        ]

# # print(
# #     m.fs.mixer_settler_ex[:]
# #     .mixer[1]
# #     .unit.mscontactor.aqueous[0.0, 1]
# #     .pH_phase["liquid"]()
# # )
# # percentage_extraction

df = pd.read_excel("data for parmest.xlsx", sheet_name="Sheet2")

element_list = ["Y", "Dy", "Gd", "Sm", "Nd", "Ce","La","Pr"]
colors = {"Y": "r", "Dy": "g", "Gd": "r", "Sm": "g", "Nd": "r", "Ce": "r","La": "g", "Pr": "g"}
pH_list = [
    m.fs.mixer_settler_ex[s]
    .mixer[1]
    .unit.mscontactor.aqueous[0.0, 1]
    .pH_phase["liquid"]()
    for s in m.system
]

# for e in element_list:
#     plt.scatter(df["pH"], df[e], marker="o", color=colors[e], label=f"{e}_exp")
#     plt.plot(
#         pH_list,
#         percentage_extraction[e],
#         color=colors[e],
#         label=f"{e}_model",
#     )
# plt.xlabel("pH")
# plt.ylabel("Extraction %")
# plt.title(f"Extraction % comparison for 5% DEHPA 10% TBP")
# plt.legend()

# plt.figure(dpi=300)
# for e in ["Y", "Dy", "La","Sm"]:
#     for i in df.index:
#         if df.loc[i,f'w_{e}'] == 1:
#             plt.scatter(df.loc[i,"pH"], df.loc[i,e], marker="o", color=colors[e], label=f"{e}_exp")
#     plt.plot(
#         pH_list,
#         percentage_extraction[e],
#         color=colors[e],
#         label=f"{e}_model",
#     )
# plt.xlabel("pH")
# plt.ylabel("Extraction %")
# plt.title(f"Extraction % comparison for 5% DEHPA 10% TBP")
# plt.legend()

fig, ax = plt.subplots(2,2,dpi=300,figsize=(7,5))
plt.suptitle('Extraction profile at 5% DEHPA 10% TBP')
for e in ["Dy", "Y"]:
    pH_exp = []
    ext = []
    for i in df.index:
        if df.loc[i, f"w_{e}"] == 1 and 4 < df.loc[i, 'dosage'] < 5:
            pH_exp.append(df.loc[i, "pH"])
            ext.append(10**df.loc[i, f'logD {e}']/(1+10**df.loc[i, f'logD {e}'])*100)
    ax[0,0].scatter(pH_exp, ext, marker="o", color=colors[e], label=f"{e}_exp")
    ax[0,0].plot(
        pH_list,
        percentage_extraction[e],
        color=colors[e],
        label=f"{e}_model",
    )
    ax[0,0].set_xlabel("pH")
    ax[0,0].set_ylabel("Extraction %")
    ax[0,0].legend(fontsize=8)
for e in ["Sm", "Gd"]:
    pH_exp = []
    ext = []
    for i in df.index:
        if df.loc[i, f"w_{e}"] == 1 and 4 < df.loc[i, 'dosage'] < 5:
            pH_exp.append(df.loc[i, "pH"])
            ext.append(10**df.loc[i, f'logD {e}']/(1+10**df.loc[i, f'logD {e}'])*100)
    ax[0,1].scatter(pH_exp, ext, marker="o", color=colors[e], label=f"{e}_exp")
    ax[0,1].plot(
        pH_list,
        percentage_extraction[e],
        color=colors[e],
        label=f"{e}_model",
    )
    ax[0,1].set_xlabel("pH")
    ax[0,1].set_ylabel("Extraction %")
    ax[0,1].legend(fontsize=8)
for e in ["Pr", "Ce"]:
    pH_exp = []
    ext = []
    for i in df.index:
        if df.loc[i, f"w_{e}"] == 1 and 4 < df.loc[i, 'dosage'] < 5:
            pH_exp.append(df.loc[i, "pH"])
            ext.append(10**df.loc[i, f'logD {e}']/(1+10**df.loc[i, f'logD {e}'])*100)
    ax[1,0].scatter(pH_exp, ext, marker="o", color=colors[e], label=f"{e}_exp")
    ax[1,0].plot(
        pH_list,
        percentage_extraction[e],
        color=colors[e],
        label=f"{e}_model",
    )
    ax[1,0].set_xlabel("pH")
    ax[1,0].set_ylabel("Extraction %")
    ax[1,0].legend(fontsize=8)
for e in ["La", "Nd"]:
    pH_exp = []
    ext = []
    for i in df.index:
        if df.loc[i, f"w_{e}"] == 1 and 4 < df.loc[i, 'dosage'] < 5:
            pH_exp.append(df.loc[i, "pH"])
            ext.append(10**df.loc[i, f'logD {e}']/(1+10**df.loc[i, f'logD {e}'])*100)
    ax[1,1].scatter(pH_exp, ext, marker="o", color=colors[e], label=f"{e}_exp")
    ax[1,1].plot(
        pH_list,
        percentage_extraction[e],
        color=colors[e],
        label=f"{e}_model",
    )
    ax[1,1].set_xlabel("pH")
    ax[1,1].set_ylabel("Extraction %")
    ax[1,1].legend(fontsize=8)
plt.tight_layout()