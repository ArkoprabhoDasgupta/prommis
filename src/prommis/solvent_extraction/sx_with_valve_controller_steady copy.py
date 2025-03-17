from pyomo.environ import (
    ConcreteModel,
    units,
)

import numpy as np

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)

from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.solvers import get_solver


from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution_new import REESolExOgParameters
from prommis.solvent_extraction.hybrid_solvent_extraction_trial import (
    SolventExtraction,
    SolventExtractionInitializer,
)


from prommis.solvent_extraction.sx_reaction_pkg import SolventExtractionReactions

m = ConcreteModel()

time_duration = 24

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()

number_of_stages = 3
stage_number = np.arange(1, number_of_stages + 1)
dosage = 5

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

m.fs.solex.mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.solex.mscontactor.volume_frac_stream[:, :, "aqueous"].fix(0.5)
m.fs.solex.cross_sec_area[:] = 1
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

# m.H_generation_term = Var(m.fs.time, stage_number)


# @m.Constraint(m.fs.time, stage_number)
# def H_generation_rule(m, t, s):
#     return m.H_generation_term[t, s] == -3 * sum(
#         m.fs.solex.mscontactor.material_transfer_term[t, s, e]
#         for e in m.fs.solex.mscontactor.stream_component_interactions
#     )


# m.fs.solex.mscontactor.aqueous_inherent_reaction_constraint[
#     :, :, "liquid", "H"
# ].deactivate()


# @m.Constraint(m.fs.time, stage_number)
# def H_reaction_rule(m, t, s):
#     return (
#         m.fs.solex.mscontactor.aqueous_inherent_reaction_generation[t, s, "liquid", "H"]
#         == m.fs.solex.mscontactor.aqueous_inherent_reaction_extent[t, s, "Ka2"]
#         + m.H_generation_term[t, s]
#     )


# m.fs.valve.control_volume.properties_out[:].pressure.fix(101235 * units.Pa)
# m.fs.valve.valve_opening[:].fix(0.5)


# m.fs.solex.mscontactor.aqueous[:, :].h2o_concentration.deactivate()

# m.fs.solex.mscontactor.material_transfer_term[:, :, "aqueous", "organic", :].fix(0)

print(dof(m))

initializer = SolventExtractionInitializer()
initializer.initialize(m.fs.solex)

# initialize_by_time_element(m.fs, m.fs.time)
solver = get_solver(solver="ipopt_v2")
solver.options["max_iter"] = 5000
solver.options["halt_on_ampl_error"] = "yes"
solver.solve(m, tee=True)

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
