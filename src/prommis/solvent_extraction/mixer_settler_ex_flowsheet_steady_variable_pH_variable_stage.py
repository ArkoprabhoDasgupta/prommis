#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import ConcreteModel, units, Set, RangeSet

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)
from idaes.core.util import to_json


from idaes.core.solvers import get_solver

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution_new import REESolExOgParameters
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package_new_modified import (
    SolventExtractionReactions,
)


m = ConcreteModel()

m.stage_number = Set(initialize=[1, 2, 3, 4, 5, 6, 7])
# m.pH_set = Set(initialize=[1, 2, 3, 4])
# pH_value = {1: 1.2, 2: 1.4, 3: 1.6, 4: 1.8, 5: 2}
dosage_value = dict(zip(m.stage_number, RangeSet(2, 8)))

m.fs = FlowsheetBlock(m.stage_number, dynamic=False)


for s in m.stage_number:

    m.fs[s].prop_o = REESolExOgParameters()
    m.fs[s].leach_soln = LeachSolutionParameters()
    m.fs[s].reaxn = SolventExtractionReactions()

    m.fs[s].mixer_settler_ex = MixerSettlerExtraction(
        number_of_stages=3,
        aqueous_stream={
            "property_package": m.fs[s].leach_soln,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs[s].prop_o,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        heterogeneous_reaction_package=m.fs[s].reaxn,
        has_holdup=True,
        settler_transformation_method="dae.finite_difference",
        settler_transformation_scheme="BACKWARD",
        settler_finite_elements=4,
    )

# dosage = 5

for s in m.stage_number:

    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "H"].fix((10**-1.5) * 1e3)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(
        (10**-1.5) * 48e3
    )
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-4)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Al"].fix(137.27)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(25.78)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(138.27)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.277)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.346)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "La"].fix(2.09)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(5)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.73)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(2.1)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.236)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.56)
    m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.09)

    m.fs[s].mixer_settler_ex.aqueous_inlet.flow_vol.fix(62.01)

    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
        975.8e3 * dosage_value[s] / 100
    )
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Al_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Ca_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Sc_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Y_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "La_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Ce_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Pr_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Nd_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Sm_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Gd_o"].fix(1e-7)
    m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Dy_o"].fix(1e-7)

    m.fs[s].mixer_settler_ex.organic_inlet.flow_vol.fix(62.01)

    m.fs[s].mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
        305.15 * units.K
    )
    m.fs[s].mixer_settler_ex.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
        305.15 * units.K
    )

    m.fs[s].mixer_settler_ex.mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)

    m.fs[s].mixer_settler_ex.organic_settler[:].unit.area.fix(1)
    m.fs[s].mixer_settler_ex.aqueous_settler[:].unit.area.fix(1)
    m.fs[s].mixer_settler_ex.aqueous_settler[:].unit.length.fix(1)
    m.fs[s].mixer_settler_ex.organic_settler[:].unit.length.fix(1)

    for e in ["Al", "Ca", "Fe", "Sc"]:
        m.fs[s].mixer_settler_ex.mixer[
            :
        ].unit.mscontactor.heterogeneous_reaction_extent[
            0.0, 1, f"{e}_mass_transfer"
        ].fix(
            0
        )
        m.fs[s].mixer_settler_ex.mixer[:].unit.distribution_extent_constraint[
            0, 1, e
        ].deactivate()

    m.fs[s].mixer_settler_ex.default_initializer().initialize(m.fs[s].mixer_settler_ex)


solver = get_solver("ipopt_v2")
results = solver.solve(m, tee=True)

percentage_extraction = {}

for e in ["La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy", "Y"]:
    percentage_extraction[e] = {}
    for s in m.stage_number:
        percentage_extraction[e][s] = (
            (
                m.fs[s].mixer_settler_ex.organic_outlet.conc_mass_comp[0, f"{e}_o"]()
                * m.fs[s].mixer_settler_ex.organic_outlet.flow_vol[0]()
                - m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, f"{e}_o"]()
                * m.fs[s].mixer_settler_ex.organic_inlet.flow_vol[0]()
            )
            / (
                m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, e]()
                * m.fs[s].mixer_settler_ex.aqueous_inlet.flow_vol[0]()
            )
        ) * 100


percentage_extraction["TREE"] = {}

for s in m.stage_number:
    percentage_extraction["TREE"][s] = (
        (
            sum(
                m.fs[s].mixer_settler_ex.organic_outlet.conc_mass_comp[0, f"{e}_o"]()
                for e in m.fs[s].leach_soln.component_list
                if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
            )
            * m.fs[s].mixer_settler_ex.organic_outlet.flow_vol[0]()
            - sum(
                m.fs[s].mixer_settler_ex.organic_inlet.conc_mass_comp[0, f"{e}_o"]()
                for e in m.fs[s].leach_soln.component_list
                if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
            )
            * m.fs[s].mixer_settler_ex.organic_inlet.flow_vol[0]()
        )
        / (
            sum(
                m.fs[s].mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, e]()
                for e in m.fs[s].leach_soln.component_list
                if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
            )
            * m.fs[s].mixer_settler_ex.aqueous_inlet.flow_vol[0]()
        )
    ) * 100


from pyomo.environ import RangeSet
import matplotlib.pyplot as plt
import numpy as np

# fig, ax = plt.subplots(3, 3, figsize=(8, 8))

# element_mapper = dict(
#     zip(RangeSet(1, 9), ["La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy", "Y", "TREE"])
# )

# for i in RangeSet(0, 2):
#     for j in RangeSet(0, 2):
#         ax[i, j].plot(
#             m.stage_number,
#             percentage_extraction[element_mapper[i + 1 + j * 3]][1],
#             label=f"{pH_value[1]}",
#             marker="o",
#         )
#         ax[i, j].plot(
#             m.stage_number,
#             percentage_extraction[element_mapper[i + 1 + j * 3]][2],
#             label=f"{pH_value[2]}",
#             marker="o",
#         )
#         ax[i, j].plot(
#             m.stage_number,
#             percentage_extraction[element_mapper[i + 1 + j * 3]][3],
#             label=f"{pH_value[3]}",
#             marker="o",
#         )
#         ax[i, j].plot(
#             m.stage_number,
#             percentage_extraction[element_mapper[i + 1 + j * 3]][4],
#             label=f"{pH_value[4]}",
#             marker="o",
#         )
#         ax[i, j].set(xlabel="Number of stages", ylabel="percent extraction")
#         ax[i, j].set_title(f"{element_mapper[i + 1 + j * 3]} recovery wrt stages")
# plt.legend(fontsize=8)
# plt.tight_layout()

categories = ["Gd", "Sm", "Ce", "TREE"]
group1 = [percentage_extraction[c][1] for c in categories]
group2 = [percentage_extraction[c][2] for c in categories]
group3 = [percentage_extraction[c][3] for c in categories]
group4 = [percentage_extraction[c][4] for c in categories]
group5 = [percentage_extraction[c][5] for c in categories]
group6 = [percentage_extraction[c][6] for c in categories]
group7 = [percentage_extraction[c][7] for c in categories]

# Set up bar positions
x = np.arange(len(categories))
width = 0.12  # Width of each bar

# Create figure and axis
fig, ax = plt.subplots(dpi=300)

# Create bars for each group
bars1 = ax.bar(x - 3 * width, group1, width, label="2 %")
bars2 = ax.bar(x - 2 * width, group2, width, label="3 %")
bars3 = ax.bar(x - 1 * width, group3, width, label="4 %")
bars4 = ax.bar(x, group4, width, label="5 %")
bars5 = ax.bar(x + 1 * width, group5, width, label="6 %")
bars6 = ax.bar(x + 2 * width, group6, width, label="7 %")
bars7 = ax.bar(x + 3 * width, group7, width, label="8 %")

# Add labels and title
ax.set_xlabel("Elements")
ax.set_ylabel("Extraction percentage (%)")
ax.set_title("Extraction profile with varying extractant dosage")
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.legend(loc="upper right", ncol=2)
ax.grid(axis="y", alpha=0.3)
plt.grid(axis="y", which="both", linewidth="0.5")
plt.minorticks_on()
plt.gca().set_axisbelow(True)

plt.tight_layout()
plt.show()
