#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import ConcreteModel, units, Set

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

m.stage_number = Set(initialize=[1, 2, 3, 4, 5, 6, 7, 8])
m.pH_set = Set(initialize=[1, 2, 3, 4])
pH_value = {1: 1.2, 2: 1.4, 3: 1.6, 4: 1.8}

m.fs = FlowsheetBlock(m.stage_number, dynamic=False)


for s in m.stage_number:

    m.fs[s].prop_o = REESolExOgParameters()
    m.fs[s].leach_soln = LeachSolutionParameters()
    m.fs[s].reaxn = SolventExtractionReactions()

    m.fs[s].mixer_settler_ex = MixerSettlerExtraction(
        m.pH_set,
        number_of_stages=s,
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

dosage = 4

for s in m.stage_number:
    for k in m.pH_set:

        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "H"].fix(
            (10 ** -pH_value[k]) * 1e3
        )
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "SO4"].fix(
            (10 ** -pH_value[k]) * 96e3
        )
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-4)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Al"].fix(422.375)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Ca"].fix(109.542)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Fe"].fix(688.266)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.032)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.124)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "La"].fix(0.986)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Ce"].fix(2.277)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.303)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Nd"].fix(0.946)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.097)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.2584)
        m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.047)

        m.fs[s].mixer_settler_ex[k].aqueous_inlet.flow_vol.fix(62.01)

        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Kerosene"].fix(
            820e3
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
            975.8e3 * dosage / 100
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Al_o"].fix(
            1.267e-5
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Ca_o"].fix(
            2.684e-5
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Fe_o"].fix(
            2.873e-6
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Sc_o"].fix(1.734)
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Y_o"].fix(2.179e-5)
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "La_o"].fix(
            0.000105
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Ce_o"].fix(0.00031)
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Pr_o"].fix(
            3.711e-5
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Nd_o"].fix(
            0.000165
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Sm_o"].fix(
            1.701e-5
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Gd_o"].fix(
            3.357e-5
        )
        m.fs[s].mixer_settler_ex[k].organic_inlet.conc_mass_comp[0, "Dy_o"].fix(
            8.008e-6
        )

        m.fs[s].mixer_settler_ex[k].organic_inlet.flow_vol.fix(62.01)

        m.fs[s].mixer_settler_ex[k].mixer[:].unit.mscontactor.aqueous[
            :, :
        ].temperature.fix(305.15 * units.K)
        m.fs[s].mixer_settler_ex[k].mixer[:].unit.mscontactor.organic[
            :, :
        ].temperature.fix(305.15 * units.K)

        m.fs[s].mixer_settler_ex[k].mixer[:].unit.mscontactor.volume[:].fix(
            0.4 * units.m**3
        )

        m.fs[s].mixer_settler_ex[k].organic_settler[:].unit.area.fix(1)
        m.fs[s].mixer_settler_ex[k].aqueous_settler[:].unit.area.fix(1)
        m.fs[s].mixer_settler_ex[k].aqueous_settler[:].unit.length.fix(1)
        m.fs[s].mixer_settler_ex[k].organic_settler[:].unit.length.fix(1)

        m.fs[s].mixer_settler_ex[k].default_initializer().initialize(
            m.fs[s].mixer_settler_ex[k]
        )

solver = get_solver("ipopt_v2")
results = solver.solve(m, tee=True)

percentage_extraction = {}

for e in ["La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy", "Y"]:
    percentage_extraction[e] = {}
    for k in m.pH_set:
        percentage_extraction[e][k] = [
            (
                (
                    m.fs[s]
                    .mixer_settler_ex[k]
                    .organic_outlet.conc_mass_comp[0, f"{e}_o"]()
                    * m.fs[s].mixer_settler_ex[k].organic_outlet.flow_vol[0]()
                    - m.fs[s]
                    .mixer_settler_ex[k]
                    .organic_inlet.conc_mass_comp[0, f"{e}_o"]()
                    * m.fs[s].mixer_settler_ex[k].organic_inlet.flow_vol[0]()
                )
                / (
                    m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, e]()
                    * m.fs[s].mixer_settler_ex[k].aqueous_inlet.flow_vol[0]()
                )
            )
            * 100
            for s in m.stage_number
        ]

percentage_extraction["TREE"] = {}

for k in m.pH_set:
    percentage_extraction["TREE"][k] = [
        (
            (
                sum(
                    m.fs[s]
                    .mixer_settler_ex[k]
                    .organic_outlet.conc_mass_comp[0, f"{e}_o"]()
                    for e in m.fs[s].leach_soln.component_list
                    if e
                    not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
                )
                * m.fs[s].mixer_settler_ex[k].organic_outlet.flow_vol[0]()
                - sum(
                    m.fs[s]
                    .mixer_settler_ex[k]
                    .organic_inlet.conc_mass_comp[0, f"{e}_o"]()
                    for e in m.fs[s].leach_soln.component_list
                    if e
                    not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
                )
                * m.fs[s].mixer_settler_ex[k].organic_inlet.flow_vol[0]()
            )
            / (
                sum(
                    m.fs[s].mixer_settler_ex[k].aqueous_inlet.conc_mass_comp[0, e]()
                    for e in m.fs[s].leach_soln.component_list
                    if e
                    not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
                )
                * m.fs[s].mixer_settler_ex[k].aqueous_inlet.flow_vol[0]()
            )
        )
        * 100
        for s in m.stage_number
    ]

from pyomo.environ import RangeSet
import matplotlib.pyplot as plt

fig, ax = plt.subplots(3, 3, figsize=(8, 8))

element_mapper = dict(
    zip(RangeSet(1, 9), ["La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy", "Y", "TREE"])
)

for i in RangeSet(0, 2):
    for j in RangeSet(0, 2):
        ax[i, j].plot(
            m.stage_number,
            percentage_extraction[element_mapper[i + 1 + j * 3]][1],
            label=f"{pH_value[1]}",
            marker="o",
        )
        ax[i, j].plot(
            m.stage_number,
            percentage_extraction[element_mapper[i + 1 + j * 3]][2],
            label=f"{pH_value[2]}",
            marker="o",
        )
        ax[i, j].plot(
            m.stage_number,
            percentage_extraction[element_mapper[i + 1 + j * 3]][3],
            label=f"{pH_value[3]}",
            marker="o",
        )
        ax[i, j].plot(
            m.stage_number,
            percentage_extraction[element_mapper[i + 1 + j * 3]][4],
            label=f"{pH_value[4]}",
            marker="o",
        )
        ax[i, j].set(xlabel="Number of stages", ylabel="percent extraction")
        ax[i, j].set_title(f"{element_mapper[i + 1 + j * 3]} recovery wrt stages")
plt.legend()
plt.tight_layout()
