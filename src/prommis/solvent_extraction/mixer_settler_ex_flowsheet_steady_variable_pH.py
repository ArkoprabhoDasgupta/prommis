#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import ConcreteModel, units

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
    m.fs.reaxn = SolventExtractionReactions()

    m.fs.reaxn.extractant_dosage = dosage

    m.fs.mixer_settler_ex = MixerSettlerExtraction(
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
    pH = 1.5
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "H"].fix((10**-pH) * 1e3)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "SO4"].fix((10**-pH) * 96e3)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-4)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Al"].fix(422.375)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(109.542)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(688.266)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.032)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.124)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "La"].fix(0.986)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(2.277)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.303)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(0.946)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.097)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.2584)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.047)

    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Al"].fix(137.27)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(25.78)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(138.27)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.277)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "La"].fix(2.09)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(5)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.73)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(2.10)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.236)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.56)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.09)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.346)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "H"].fix(
        10 ** (-pH) * units.gram / units.L
    )
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(
        10 ** (-pH) * 48 * units.gram / units.L
    )
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-4)

    m.fs.mixer_settler_ex.aqueous_inlet.flow_vol.fix(62.01)

    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
        975.8e3 * dosage / 100
    )
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Al_o"].fix(1.267e-5)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Ca_o"].fix(2.684e-5)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(2.873e-6)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Sc_o"].fix(1.734)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Y_o"].fix(2.179e-5)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "La_o"].fix(0.000105)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Ce_o"].fix(0.00031)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Pr_o"].fix(3.711e-5)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Nd_o"].fix(0.000165)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Sm_o"].fix(1.701e-5)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Gd_o"].fix(3.357e-5)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Dy_o"].fix(8.008e-6)

    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Al_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Ca_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Sc_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Y_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "La_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Ce_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Pr_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Nd_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Sm_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Gd_o"].fix(1e-7)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, "Dy_o"].fix(1e-7)

    m.fs.mixer_settler_ex.organic_inlet.flow_vol.fix(62.01)

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
        305.15 * units.K
    )
    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
        305.15 * units.K
    )

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)

    m.fs.mixer_settler_ex.organic_settler[:].unit.area.fix(1)
    m.fs.mixer_settler_ex.aqueous_settler[:].unit.area.fix(1)
    m.fs.mixer_settler_ex.aqueous_settler[:].unit.length.fix(1)
    m.fs.mixer_settler_ex.organic_settler[:].unit.length.fix(1)

    for e in ["Al", "Ca", "Fe", "Sc"]:
        m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.heterogeneous_reaction_extent[
            0.0, 1, f"{e}_mass_transfer"
        ].fix(0)
        m.fs.mixer_settler_ex.mixer[:].unit.distribution_extent_constraint[
            0, 1, e
        ].deactivate()


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
number_of_stages = 3

if __name__ == "__main__":
    m, results = main(dosage, number_of_stages)

percentage_extraction = {}

# for s in m.scenario:
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]:
        percentage_extraction[e] = [
            (
                (
                    m.fs.mixer_settler_ex.organic_outlet.conc_mass_comp[0, f"{e}_o"]()
                    * m.fs.mixer_settler_ex.organic_outlet.flow_vol[0]()
                    - m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, f"{e}_o"]()
                    * m.fs.mixer_settler_ex.organic_inlet.flow_vol[0]()
                )
                / (
                    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, e]()
                    * m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[0]()
                )
            )
            * 100
        ]

percentage_extraction["TREE"] = [
    (
        (
            sum(
                m.fs.mixer_settler_ex.organic_outlet.conc_mass_comp[0, f"{e}_o"]()
                for e in m.fs.leach_soln.component_list
                if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
            )
            * m.fs.mixer_settler_ex.organic_outlet.flow_vol[0]()
            - sum(
                m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[0, f"{e}_o"]()
                for e in m.fs.leach_soln.component_list
                if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
            )
            * m.fs.mixer_settler_ex.organic_inlet.flow_vol[0]()
        )
        / (
            sum(
                m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, e]()
                for e in m.fs.leach_soln.component_list
                if e not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
            )
            * m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[0]()
        )
    )
    * 100
]

bar_x = []
bar_y = []
for k, v in percentage_extraction.items():
    bar_x.append(k)
    bar_y.append(v[0])

import matplotlib.pyplot as plt

plt.figure(dpi=300)
plt.bar(bar_x, bar_y)
plt.xlabel("Elements")
plt.ylabel("Percentage extraction (%)")
plt.title("Extraction profile, stages = 3, feed pH = 1.5, DEHPA dosage = 5%")
plt.grid(axis="y", which="both", linewidth="0.5")
plt.minorticks_on()
plt.gca().set_axisbelow(True)
