from pyomo.environ import ConcreteModel, units, TransformationFactory, Var
from pyomo.dae.flatten import flatten_dae_components

import numpy as np

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)

from idaes.core.solvers import get_solver
from idaes.core.util import from_json

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)
from idaes.core.util import DiagnosticsToolbox
import matplotlib.pyplot as plt


def build_model(time_duration, dosage, number_of_stages):
    """
    Build model with dosage and number of stages as arguments
    """
    m = ConcreteModel()

    m.fs = FlowsheetBlock(
        dynamic=True, time_set=[0, time_duration], time_units=units.hour
    )

    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.reaxn = SolventExtractionReactions()

    m.fs.reaxn.extractant_dosage = dosage

    m.fs.solex = SolventExtraction(
        number_of_finite_elements=number_of_stages,
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
        reaction_package=m.fs.reaxn,
        has_holdup=True,
    )

    return m


def discretization(m):
    """
    Discretization of the time domain
    """

    m.discretizer = TransformationFactory("dae.collocation")
    m.discretizer.apply_to(m, nfe=10, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")

    # m.discretizer = TransformationFactory("dae.finite_difference")
    # m.discretizer.apply_to(m, nfe=7, wrt=m.fs.time, scheme="BACKWARD")


def copy_first_steady_state(m):
    """
    Function that propagates initial steady state guess to future time points.
    This function is used to initialize all the time discrete variables to the
    initial steady state value.
    """
    regular_vars, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
    # Copy initial conditions forward
    for var in time_vars:
        for t in m.fs.time:
            if t == m.fs.time.first():
                continue
            else:
                var[t].value = var[m.fs.time.first()].value


def set_inputs(m, dosage):
    """
    Set feed values to the model
    """

    m.fs.solex.mscontactor.volume[:].fix(0.4 * units.m**3)
    m.fs.solex.mscontactor.volume_frac_stream[:, :, "aqueous"].fix(0.5)
    m.fs.solex.area_cross_stage[:] = 1
    m.fs.solex.elevation[:] = 0

    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H2O"].fix(1e6)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H"].fix(10.75)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["SO4"].fix(100)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["HSO4"].fix(1e4)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Al"].fix(422.375)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ca"].fix(109.542)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Cl"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Fe"].fix(688.266)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sc"].fix(0.032)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Y"].fix(0.124)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["La"].fix(0.986)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ce"].fix(2.277)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Pr"].fix(0.303)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Nd"].fix(0.946)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sm"].fix(0.097)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Gd"].fix(0.2584)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Dy"].fix(0.047)

    m.fs.solex.mscontactor.aqueous_inlet_state[:].flow_vol.fix(62.01)

    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Kerosene"].fix(820e3)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["DEHPA"].fix(
        975.8e3 * dosage / 100
    )
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Al_o"].fix(1.267e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ca_o"].fix(2.684e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Fe_o"].fix(2.873e-6)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sc_o"].fix(1.734)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Y_o"].fix(2.179e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["La_o"].fix(0.000105)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ce_o"].fix(0.00031)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Pr_o"].fix(3.711e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Nd_o"].fix(0.000165)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sm_o"].fix(1.701e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Gd_o"].fix(3.357e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Dy_o"].fix(8.008e-6)

    m.fs.solex.mscontactor.organic_inlet_state[:].flow_vol.fix(62.01)

    m.fs.solex.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].temperature.fix(305.15 * units.K)
    m.fs.solex.mscontactor.organic[:, :].temperature.fix(305.15 * units.K)
    m.fs.solex.mscontactor.organic_inlet_state[:].temperature.fix(305.15 * units.K)

    for e in m.fs.leach_soln.component_list:
        if e not in ["H2O", "HSO4"]:
            m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp[e].fix()

    m.fs.solex.mscontactor.aqueous[0, :].flow_vol.fix()

    m.fs.solex.mscontactor.aqueous_inherent_reaction_extent[0.0, :, "Ka2"].fix()

    m.fs.solex.mscontactor.organic[0, :].flow_vol.fix()
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["DEHPA"].fix()

    for e in m.fs.reaxn.element_list:
        m.fs.solex.mscontactor.heterogeneous_reaction_extent[
            0.0, :, f"{e}_mass_transfer"
        ].fix()


def main(dosage, number_of_stages, time_duration):
    """
    Main function to run the model
    """

    m = build_model(time_duration, dosage, number_of_stages)
    discretization(m)
    from_json(m, fname="solvent_extraction.json")
    copy_first_steady_state(m)
    set_inputs(m, dosage)

    return m


"""
Build the model, initialize and solve it.
"""

dosage = 5
number_of_stages = 3
stage_number = np.arange(1, number_of_stages + 1)

time_duration = 12

m = main(dosage, number_of_stages, time_duration)

perturb_time = 4
for t in m.fs.time:
    if t <= perturb_time:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].flow_vol.fix(62.01)
    else:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].flow_vol.fix(72.01)

solver = get_solver(solver="ipopt_v2")
solver.solve(m, tee=True)
percentage_recovery = {}
REE_set = ["Y", "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy"]

for e in REE_set:
    percentage_recovery[e] = [
        (
            1
            - (
                m.fs.solex.mscontactor.aqueous[t, number_of_stages].conc_mass_comp[e]()
                * m.fs.solex.mscontactor.aqueous[t, number_of_stages].flow_vol()
            )
            / (
                m.fs.solex.mscontactor.aqueous_inlet_state[t].conc_mass_comp[e]()
                * m.fs.solex.mscontactor.aqueous_inlet_state[t].flow_vol()
            )
        )
        * 100
        for t in m.fs.time
    ]


for e in REE_set:
    plt.plot(m.fs.time, m.fs.solex.mscontactor.organic[:, 1].conc_mass_comp[f"{e}_o"]())
plt.legend(REE_set)
plt.axvline(x=perturb_time, color="red", linestyle="--")
