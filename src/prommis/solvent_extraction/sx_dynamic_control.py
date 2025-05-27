from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    units,
    TransformationFactory,
    Var,
    value,
    log10,
    Suffix,
    Constraint,
)
from pyomo.dae.flatten import flatten_dae_components
from pyomo.network import Arc

import numpy as np

import matplotlib.pyplot as plt

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
import idaes.core.solvers.petsc as petsc
from idaes.core.util.initialization import initialize_by_time_element
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util import from_json, DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.scaling import (
    set_scaling_factor,
    CustomScalerBase,
    report_scaling_factors,
)
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Valve
from idaes.models.control.controller import (
    PIDController,
    ControllerType,
    ControllerMVBoundType,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)

m = ConcreteModel()

time_duration = 24

m.fs = FlowsheetBlock(dynamic=True, time_set=[0, time_duration], time_units=units.hour)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()

dosage = 5
number_of_stages = 1
stage_number = np.arange(1, number_of_stages + 1)

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


def _valve_pressure_flow_cb(b):

    b.Cv = Var(initialize=1e-2, units=units.dm**2)
    b.Cv.fix()

    @b.Constraint(b.flowsheet().time)
    def pressure_flow_equation(b, t):
        rho_aqueous = units.convert(
            b.control_volume.properties_in[t].dens_mass,
            to_units=units.kg / (units.m**3),
        )
        # rho_aqueous = units.convert(
        #     sum(
        #         b.control_volume.properties_out[t].conc_mass_comp[p]
        #         for p in m.fs.leach_soln.component_list
        #     ),
        #     to_units=units.kg / (units.m**3),
        # )
        Po = units.convert(
            b.control_volume.properties_out[t].pressure, to_units=units.Pa
        )
        Pi = units.convert(
            b.control_volume.properties_in[t].pressure, to_units=units.Pa
        )
        vel_head = units.convert(
            (Pi - Po) / rho_aqueous, to_units=(units.dm**2 / units.hr**2)
        )

        F = units.convert(
            b.control_volume.properties_out[t].flow_vol,
            to_units=units.L / units.hr,
        )
        Cv = b.Cv
        fun = b.valve_function[t]
        return F**2 == ((Cv * fun) ** 2) * vel_head
        # return (Cv * F) ** 2 == vel_head * fun**2


m.fs.valve = Valve(
    dynamic=False,
    has_holdup=False,
    material_balance_type=MaterialBalanceType.componentTotal,
    property_package=m.fs.leach_soln,
    pressure_flow_callback=_valve_pressure_flow_cb,
)


m.fs.sx_to_v = Arc(
    source=m.fs.solex.mscontactor.aqueous_outlet, destination=m.fs.valve.inlet
)

m.fs.control = PIDController(
    process_var=m.fs.solex.mscontactor.volume_frac_stream[:, 1, "aqueous"],
    manipulated_var=m.fs.valve.valve_opening,
    controller_type=ControllerType.PI,
)

TransformationFactory("network.expand_arcs").apply_to(m.fs)

m.discretizer = TransformationFactory("dae.finite_difference")
m.discretizer.apply_to(m, nfe=40, wrt=m.fs.time, scheme="BACKWARD")

from_json(m, fname="solvent_extraction.json")


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


copy_first_steady_state(m)

# Fixing inlet and feed stuffs

m.fs.solex.mscontactor.volume[:].fix(100 * units.L)
# for t in m.fs.time:
#     m.fs.solex.mscontactor.volume_frac_stream[t, :, "aqueous"].fix(0.5 + 0.3 * (t / 24))

m.fs.solex.mscontactor.volume_frac_stream[0, :, "aqueous"].fix(0.3)
m.fs.solex.mscontactor.volume_frac_stream.setub(1)
m.fs.solex.mscontactor.volume_frac_stream.setlb(0)

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


m.fs.solex.mscontactor.aqueous_inlet_state[:].temperature.fix(305.15 * units.K)
m.fs.solex.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)

m.fs.control.gain_p.fix(3)
m.fs.control.gain_i.fix(1)
for t in m.fs.time:
    if 0 <= t <= 4:
        m.fs.control.setpoint[t].fix(0.3)
    elif 4 < t <= 15:
        m.fs.control.setpoint[t].fix(0.5)
    else:
        m.fs.control.setpoint[t].fix(0.6)
# m.fs.control.setpoint.fix(0.5)
m.fs.control.mv_ref.fix(0)
# m.fs.control.derivative_term[0].fix(1e-4)
# m.fs.control.mv_eqn[:].deactivate()

m.fs.valve.control_volume.properties_out[:].pressure.fix(101235 * units.Pa)
m.fs.valve.control_volume.properties_out[:].temperature.fix(305.15 * units.K)
m.fs.valve.control_volume.enthalpy_balances[:].deactivate()
m.fs.valve.valve_opening[:].unfix()

print(dof(m))

solver = get_solver(solver="ipopt_v2")
solver.solve(m, tee=True)
