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
from idaes.core.util.misc import add_object_reference
import idaes.core.solvers.petsc as petsc
from idaes.core.util.initialization import initialize_by_time_element
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util import from_json, DiagnosticsToolbox, StoreSpec
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
from prommis.solvent_extraction.solvent_extraction_reaction_package_new_modified import (
    SolventExtractionReactions,
)
from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution_new import REESolExOgParameters
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
)


time_duration = 24
number_of_stages = 1
dosage = 5


m = ConcreteModel()

# make the mixer settler ex
m.fs = FlowsheetBlock(
    dynamic=True, time_set=[0, time_duration], time_units=units.hour
)
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


# m.fs.mixer_settler_ex.mixer[:].unit.volume_fraction_constraint.deactivate()

# m.discretizer = TransformationFactory("dae.collocation")
# m.discretizer.apply_to(m, nfe=2, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")

def _valve_pressure_flow_cb(b):

    b.Cv = Var(initialize=1e-2, units=units.dm**2)
    b.Cv.fix()

    @b.Constraint(b.flowsheet().time)
    def pressure_flow_equation(b, t):
        # rho_aqueous = units.convert(
        #     b.control_volume.properties_in[t].dens_mass,
        #     to_units=units.kg / (units.m**3),
        # )
        rho_aqueous = units.convert(
            sum(
                b.control_volume.properties_out[t].conc_mass_comp[p]
                for p in m.fs.leach_soln.component_list
            ),
            to_units=units.kg / (units.m**3),
        )
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
    source=m.fs.mixer_settler_ex.aqueous_outlet, destination=m.fs.valve.inlet
)

m.fs.control = PIDController(
    process_var=m.fs.mixer_settler_ex.mixer[1].unit.mscontactor.volume_frac_stream[:, 1, "aqueous"],
    manipulated_var=m.fs.valve.valve_opening,
    controller_type=ControllerType.PI,
)

TransformationFactory("network.expand_arcs").apply_to(m.fs)

# m.discretizer = TransformationFactory("dae.finite_difference")
# m.discretizer.apply_to(m, nfe=6, wrt=m.fs.time, scheme="BACKWARD")

m.discretizer = TransformationFactory("dae.collocation")
m.discretizer.apply_to(m, nfe=6, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")


path_name="mixer_settler_extraction.json"
from_json(m, fname=path_name, wts=StoreSpec.value())

m.fs.mixer_settler_ex.mixer[:].unit.volume_fraction_constraint.deactivate()

regular_vars, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
# Copy initial conditions forward
for var in time_vars:
    for t in m.fs.time:
        if t == m.fs.time.first():
            continue
        else:
            var[t].value = var[m.fs.time.first()].value

m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "H2O"].fix(1e6)
# m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "H"].fix(10.75)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "SO4"].fix(100)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "HSO4"].fix(1e4)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Al"].fix(422.375)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Ca"].fix(109.542)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Cl"].fix(1e-7)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Fe"].fix(688.266)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Sc"].fix(0.032)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Y"].fix(0.124)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "La"].fix(0.986)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Ce"].fix(2.277)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Pr"].fix(0.303)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Nd"].fix(0.946)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Sm"].fix(0.097)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Gd"].fix(0.2584)
m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Dy"].fix(0.047)

perturb_time = 10
for t in m.fs.time:
    if t <= perturb_time:
        m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t].fix(62.01)
    else:
        m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t].fix(62.01)
        # m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t].fix(62.01)
    if t <= perturb_time:
        m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Gd"].fix(0.2584)
    else:
        # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Gd"].fix(0.2584)
        m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Gd"].fix(0.2584)
    if t <= perturb_time * 3:
        m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"].fix(10.75)
    # elif perturb_time <= t < perturb_time * 2:
    #     m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"].fix(8.75)
    else:
        m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"].fix(10.75)
        # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"].fix(10.75)

m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Kerosene"].fix(820e3)
# m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "DEHPA"].fix(
#     975.8e3 * dosage / 100
# )
for t in m.fs.time:
    if t <= perturb_time * 3:
        m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[t, "DEHPA"].fix(
            975.8e3 * dosage / 100
        )
    else:
        m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[t, "DEHPA"].fix(
            975.8e3 * dosage / 100
        )
        # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[t, "DEHPA"].fix(
        #     975.8e3 * dosage / 100
        # )
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Al_o"].fix(1.267e-5)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Ca_o"].fix(2.684e-5)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Fe_o"].fix(2.873e-6)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Sc_o"].fix(1.734)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Y_o"].fix(2.179e-5)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "La_o"].fix(0.000105)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Ce_o"].fix(0.00031)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Pr_o"].fix(3.711e-5)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Nd_o"].fix(0.000165)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Sm_o"].fix(1.701e-5)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Gd_o"].fix(3.357e-5)
m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Dy_o"].fix(8.008e-6)

m.fs.mixer_settler_ex.organic_inlet.flow_vol.fix(62.01)

# Fixing mixer parameters

m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)

m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
    305.15 * units.K
)

# Fixing settler parameters

m.fs.mixer_settler_ex.organic_settler[:].unit.area.fix(1)
m.fs.mixer_settler_ex.aqueous_settler[:].unit.area.fix(1)
m.fs.mixer_settler_ex.aqueous_settler[:].unit.length.fix(0.1)
m.fs.mixer_settler_ex.organic_settler[:].unit.length.fix(0.1)


for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "HSO4"]:
        m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[
            0, :
        ].conc_mass_comp[e].fix()

m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.volume_frac_stream[
    0, :, "aqueous"
].fix()
m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[0, :].flow_vol.fix()

m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous_inherent_reaction_extent[
    0.0, :, "Ka2"
].fix()

m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[0, :].flow_vol.fix()
m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[0, :].conc_mass_comp[
    "DEHPA"
].fix()

for e in m.fs.reaxn.element_list:
    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.heterogeneous_reaction_extent[
        0.0, :, f"{e}_mass_transfer"
    ].fix()

for e in ["Al", "Ca", "Fe", "Sc"]:
    for s in m.fs.mixer_settler_ex.elements:
        m.fs.mixer_settler_ex.mixer[
            s
        ].unit.mscontactor.heterogeneous_reaction_extent[
            :, :, f"{e}_mass_transfer"
        ].fix(
            0
        )
        m.fs.mixer_settler_ex.mixer[s].unit.mscontactor.organic[
            0.0, 1
        ].conc_mass_comp[f"{e}_o"].fix()
        m.fs.mixer_settler_ex.mixer[s].unit.distribution_extent_constraint[
            :, :, e
        ].deactivate()

# set variable values in the settler at t=0

for s in m.fs.mixer_settler_ex.elements:
    for x in m.fs.mixer_settler_ex.aqueous_settler[s].unit.length_domain:
        if x != 0:
            for e in m.fs.leach_soln.component_list:
                if e not in ["H2O", "HSO4"]:
                    m.fs.mixer_settler_ex.aqueous_settler[s].unit.properties[
                        0, x
                    ].conc_mass_comp[e].fix()
            m.fs.mixer_settler_ex.aqueous_settler[s].unit.properties[
                0, x
            ].flow_vol.fix()
            m.fs.mixer_settler_ex.aqueous_settler[s].unit.inherent_reaction_extent[
                0, x, "Ka2"
            ].fix()

    for x in m.fs.mixer_settler_ex.organic_settler[s].unit.length_domain:
        if x != 0:
            for e in m.fs.prop_o.component_list:
                if e not in ["Kerosene"]:
                    m.fs.mixer_settler_ex.organic_settler[s].unit.properties[
                        0, x
                    ].conc_mass_comp[e].fix()
            m.fs.mixer_settler_ex.organic_settler[s].unit.properties[
                0, x
            ].flow_vol.fix()



m.fs.control.gain_p.fix(30)
m.fs.control.gain_i.fix(10)
# m.fs.control.gain_d.fix(0)
for t in m.fs.time:
    if t <= 10:
        m.fs.control.setpoint[t].fix(0.5)
    else:
        m.fs.control.setpoint[t].fix(0.7)
# m.fs.control.setpoint.fix(0.5)
m.fs.control.mv_ref.fix(0)
# m.fs.control.derivative_term[0].fix(1e-4)
# m.fs.control.mv_eqn[:].deactivate()

m.fs.valve.control_volume.properties_out[:].pressure.fix(101235 * units.Pa)
m.fs.valve.control_volume.properties_out[:].temperature.fix(305.15 * units.K)
m.fs.valve.control_volume.enthalpy_balances[:].deactivate()
# m.fs.valve.control_volume.work[:].fix(0)
m.fs.valve.valve_opening[:].unfix()
# m.fs.valve.valve_opening[:].fix(0.8)


print(dof(m))

# assert 1==2

# initialize_by_time_element(m.fs, m.fs.time)
solver = get_solver(solver="ipopt_v2")
# solver.options["max_iter"] = 10000
# solver.solve(m, tee=True)

# report_scaling_factors(m, descend_into=True)

# solver = get_solver(
#     "ipopt_v2", writer_config={"linear_presolve": True, "scale_model": True}
# )
solver.options["max_iter"] = 10000
solver.solve(m, tee=True)
