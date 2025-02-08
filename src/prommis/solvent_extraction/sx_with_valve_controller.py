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
from prommis.solvent_extraction.hybrid_solvent_extraction import SolventExtraction
from prommis.solvent_extraction.D_variant_model import D_calculation

m = ConcreteModel()

time_duration = 24

Elements = ["Y", "Ce", "Nd", "Sm", "Gd", "Dy"]

m.fs = FlowsheetBlock(dynamic=True, time_set=[0, time_duration], time_units=units.hour)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()

number_of_stages = 3
stage_number = np.arange(1, number_of_stages + 1)

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
)


# def _valve_pressure_flow_cb(b):

#     b.Cv = Var(initialize=1e5)
#     b.Cv.fix()

#     @b.Constraint(b.flowsheet().time)
#     def pressure_flow_equation(b, t):
#         rho_aqueous = units.convert(
#             b.control_volume.properties_in[t].dens_mass,
#             to_units=units.kg / (units.m**3),
#         )
#         Po = units.convert(
#             b.control_volume.properties_out[t].pressure, to_units=units.Pa
#         )
#         Pi = units.convert(
#             b.control_volume.properties_in[t].pressure, to_units=units.Pa
#         )
#         F = units.convert(
#             b.control_volume.properties_in[t].flow_vol,
#             to_units=(units.m**3) / units.sec,
#         )
#         Cv = b.Cv
#         fun = b.valve_function[t]
#         # return F**2 == (((Cv * fun) ** 2) * (Pi - Po)) / rho_aqueous
#         return (Cv * F) ** 2 == ((fun**2) * (Pi - Po)) / rho_aqueous


# m.fs.valve = Valve(
#     dynamic=False,
#     has_holdup=False,
#     material_balance_type=MaterialBalanceType.componentTotal,
#     property_package=m.fs.leach_soln,
#     pressure_flow_callback=_valve_pressure_flow_cb,
# )


# m.fs.sx_to_v = Arc(
#     source=m.fs.solex.mscontactor.aqueous_outlet, destination=m.fs.valve.inlet
# )

# m.fs.control = PIDController(
#     process_var=m.fs.solex.mscontactor.volume_frac_stream[:, 1, "aqueous"],
#     manipulated_var=m.fs.valve.valve_opening,
#     controller_type=ControllerType.P,
# )

TransformationFactory("network.expand_arcs").apply_to(m.fs)

m.discretizer = TransformationFactory("dae.finite_difference")
m.discretizer.apply_to(m, nfe=30, wrt=m.fs.time, scheme="BACKWARD")


# Fixing inlet and feed stuffs

m.fs.solex.mscontactor.volume[:].fix(100 * units.L)
# for t in m.fs.time:
#     m.fs.solex.mscontactor.volume_frac_stream[t, :, "aqueous"].fix(0.5 + 0.3 * (t / 24))

m.fs.solex.mscontactor.volume_frac_stream[:, :, "aqueous"].fix(0.5)
m.fs.solex.mscontactor.volume_frac_stream.setub(1)
m.fs.solex.mscontactor.volume_frac_stream.setlb(0)
m.pH = Var(m.fs.time, stage_number)


@m.Constraint(m.fs.time, stage_number)
def pH_value(m, t, s):
    return m.pH[t, s] == -log10(
        m.fs.solex.mscontactor.aqueous[t, s].conc_mol_comp["H"] * units.L / units.mol
    )
    # return (
    #     10 ** (-m.pH[t, s]) * units.mol / units.L
    #     == m.fs.solex.mscontactor.aqueous[t, s].conc_mol_comp["H"]
    # )


@m.Constraint(m.fs.time, stage_number, Elements)
def distribution_calculation(m, t, s, e):
    a, b = D_calculation(e, 5)
    return m.fs.solex.distribution_coefficient[t, s, "aqueous", "organic", e] == 10 ** (
        a * m.pH[t, s] + b
    )


@m.Constraint(m.fs.time, stage_number)
def pressure_calculation(m, t, s):
    g = 9.8 * (units.m) / units.sec**2
    P_atm = 101325 * units.Pa
    A = 1 * units.m**2
    P_aq = units.convert(
        (
            m.fs.leach_soln.dens_mass
            * g
            * m.fs.solex.mscontactor.volume[s]
            * m.fs.solex.mscontactor.volume_frac_stream[t, s, "aqueous"]
            / A
        ),
        to_units=units.Pa,
    )
    P_org = units.convert(
        (
            m.fs.prop_o.dens_mass
            * g
            * m.fs.solex.mscontactor.volume[s]
            * m.fs.solex.mscontactor.volume_frac_stream[t, s, "organic"]
            / A
        ),
        to_units=units.Pa,
    )
    return m.fs.solex.mscontactor.aqueous[t, s].pressure == P_aq + P_org + P_atm


for s in stage_number:
    if s == 1:
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 5.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 3 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 24.7 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 99.1 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 32.4 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 58.2 / 100

    else:
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 4.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 12.3 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 6.4 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 16.7 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 23.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 15.1 / 100

m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H2O"].fix(1e6)
# m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H"].fix(10.75)
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
# m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Dy"].fix(0.047)

for t in m.fs.time:
    if t <= 10:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].conc_mass_comp["H"].fix(10.75)
    else:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].conc_mass_comp["H"].fix(10.75)

for t in m.fs.time:
    if t <= 10:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].conc_mass_comp["Dy"].fix(0.047)
    else:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].conc_mass_comp["Dy"].fix(0.047)


for t in m.fs.time:
    if t <= 10:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].flow_vol.fix(62.01)
    else:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].flow_vol.fix(62.01)


# m.fs.solex.mscontactor.aqueous_inlet_state[:].flow_vol.fix(62.01)

m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["DEHPA"].fix(975.8e3)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Al"].fix(1.267e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ca"].fix(2.684e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Fe"].fix(2.873e-6)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sc"].fix(1.734)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Y"].fix(2.179e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["La"].fix(0.000105)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ce"].fix(0.00031)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Pr"].fix(3.711e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Nd"].fix(0.000165)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sm"].fix(1.701e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Gd"].fix(3.357e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Dy"].fix(8.008e-6)

for t in m.fs.time:
    if t <= 10:
        m.fs.solex.mscontactor.organic_inlet_state[t].flow_vol.fix(62.01)
    else:
        m.fs.solex.mscontactor.organic_inlet_state[t].flow_vol.fix(62.01)

# m.fs.solex.mscontactor.organic_inlet_state[:].flow_vol.fix(62.01)

m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["H"].fix(10.75)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["SO4"].fix(100)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["HSO4"].fix(1e4)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Al"].fix(422.375)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ca"].fix(109.542)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Cl"].fix(1e-7)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Fe"].fix(688.266)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sc"].fix(0.32)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Y"].fix(0.124)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["La"].fix(0.986)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ce"].fix(2.277)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Pr"].fix(0.303)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Nd"].fix(0.946)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sm"].fix(0.097)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Gd"].fix(0.2584)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Dy"].fix(0.047)

m.fs.solex.mscontactor.aqueous_inherent_reaction_extent[0, :, "Ka2"].fix(1e-8)
m.fs.solex.mscontactor.aqueous[0, :].flow_vol.fix(62.01)

m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Al"].fix(1.267e-5)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ca"].fix(2.684e-5)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Fe"].fix(2.873e-6)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sc"].fix(1.734)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Y"].fix(2.179e-5)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["La"].fix(0.000105)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ce"].fix(0.00031)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Pr"].fix(3.711e-5)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Nd"].fix(0.000165)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sm"].fix(1.701e-5)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Gd"].fix(3.357e-5)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Dy"].fix(8.008e-6)

m.fs.solex.mscontactor.organic[0, :].flow_vol.fix(62.01)

m.fs.solex.mscontactor.aqueous[0.0, :].hso4_dissociation.deactivate()

for e in Elements:
    m.fs.solex.mass_transfer_constraint[0, :, "aqueous", "organic", e].deactivate()

for e in Elements:
    m.fs.solex.mscontactor.material_transfer_term[0.0, :, "aqueous", "organic", e].fix(
        -1
    )

m.fs.solex.mscontactor.aqueous_inlet_state[:].temperature.fix(305.15 * units.K)
m.fs.solex.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)

m.H_generation_term = Var(m.fs.time, stage_number)


@m.Constraint(m.fs.time, stage_number)
def H_generation_rule(m, t, s):
    return m.H_generation_term[t, s] == -3 * sum(
        m.fs.solex.mscontactor.material_transfer_term[t, s, e]
        for e in m.fs.solex.mscontactor.stream_component_interactions
    )


m.fs.solex.mscontactor.aqueous_inherent_reaction_constraint[
    :, :, "liquid", "H"
].deactivate()


@m.Constraint(m.fs.time, stage_number)
def H_reaction_rule(m, t, s):
    return (
        m.fs.solex.mscontactor.aqueous_inherent_reaction_generation[t, s, "liquid", "H"]
        == m.fs.solex.mscontactor.aqueous_inherent_reaction_extent[t, s, "Ka2"]
        + m.H_generation_term[t, s]
    )


# m.fs.control.gain_p.fix(10)
# m.fs.control.gain_i.fix(3)
# m.fs.control.gain_d.fix(0)
# for t in m.fs.time:
#     if t <= 10:
#         m.fs.control.setpoint[t].fix(0.5)
#     else:
#         m.fs.control.setpoint[t].fix(0.5)
# m.fs.control.setpoint.fix(0.5)
# m.fs.control.mv_ref.fix(0)
# m.fs.control.derivative_term[0].fix(1e-4)
# m.fs.control.mv_eqn[:].deactivate()

# m.fs.valve.control_volume.properties_out[:].pressure.fix(101235 * units.Pa)
# m.fs.valve.valve_opening[:].unfix()
# m.fs.valve.valve_opening[:].fix(0.5)

print(dof(m))

# initialize_by_time_element(m.fs, m.fs.time)
solver = get_solver(solver="ipopt_v2")
solver.options["max_iter"] = 5000
solver.solve(m, tee=True)

percent_recovery = {}
for ei, e in enumerate(Elements):
    for si, s in enumerate(stage_number):
        percent_recovery[e, s] = [
            (
                1
                - (
                    (
                        m.fs.solex.mscontactor.aqueous[t, s].conc_mass_comp[e]()
                        * m.fs.solex.mscontactor.aqueous[t, s].flow_vol()
                    )
                    / (
                        m.fs.solex.mscontactor.aqueous_inlet_state[t].conc_mass_comp[
                            e
                        ]()
                        * m.fs.solex.mscontactor.aqueous_inlet_state[t].flow_vol()
                    )
                )
            )
            * 100
            for t in m.fs.time
        ]


for e in Elements:
    plt.plot(m.fs.time, percent_recovery[e, 1])
plt.legend(Elements)
plt.xlabel("time, hrs")
plt.ylabel("percent recovery, %")
plt.title("Aqueous phase percent recovery graph w.r.t. time")

# result = petsc.petsc_dae_by_time_element(
#     m,
#     time=m.fs.time,
#     ts_options={
#         "--ts_type": "beuler",
#         "--ts_dt": 0.1,
#         "--ts_monitor": "",  # set initial step to 0.1
#         "--ts_save_trajectory": 1,
#     },
# )
# tj = result.trajectory  # trajectroy data
