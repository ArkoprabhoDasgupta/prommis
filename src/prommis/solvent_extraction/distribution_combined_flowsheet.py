from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory
from pyomo.network import Arc

import numpy as np

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from idaes.core.util.model_statistics import degrees_of_freedom as dof


m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()

# Loading model

m.fs.solex_loading = SolventExtraction(
    number_of_finite_elements=1,
    dynamic=False,
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
    aqueous_to_organic=True
)

# Stripping model

m.fs.solex_stripping = SolventExtraction(
    number_of_finite_elements=1,
    dynamic=False,
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
    aqueous_to_organic=False
)

# Distribution coefficient values for the different operations

Elements = ["Y", "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy"]
slope = dict(zip(Elements, [2.39, 0.54, 0.55, 0.29, 0.5, 0.7, 1.12, 2.03]))
intercept = dict(zip(Elements, [-2.32, -1.93, -1.94, -1.48, -1.89, -2, -2.22, -2.44]))

pH_loading = 1.524
pH_stripping = 0.241

for e in Elements:
    m.fs.solex_loading.distribution_coefficient[:, "aqueous", "organic", e] = 10**(slope[e]*pH_loading + intercept[e])
    m.fs.solex_stripping.distribution_coefficient[:, "aqueous", "organic", e] = 10**(slope[e]*pH_stripping + intercept[e])

m.fs.solex_loading.distribution_coefficient[:, "aqueous", "organic", "Sc"] = 1
m.fs.solex_loading.distribution_coefficient[:, "aqueous", "organic", "Al"] = 1
m.fs.solex_loading.distribution_coefficient[:, "aqueous", "organic", "Fe"] = 1
m.fs.solex_loading.distribution_coefficient[:, "aqueous", "organic", "Ca"] = 1
m.fs.solex_stripping.distribution_coefficient[:, "aqueous", "organic", "Sc"] = 1
m.fs.solex_stripping.distribution_coefficient[:, "aqueous", "organic", "Al"] = 1
m.fs.solex_stripping.distribution_coefficient[:, "aqueous", "organic", "Fe"] = 1
m.fs.solex_stripping.distribution_coefficient[:, "aqueous", "organic", "Ca"] = 1

# Connecting loading outlet and stripping inlet

m.loading_to_stripping = Arc(source = m.fs.solex_loading.mscontactor.organic_outlet, destination=m.fs.solex_stripping.mscontactor.organic_inlet)
TransformationFactory("network.expand_arcs").apply_to(m)

# Fixing feed connections

eps = 1e-18

m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H2O"].fix(1e6)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H"].fix(eps)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["SO4"].fix(eps)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["HSO4"].fix(eps)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Al"].fix(eps)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ca"].fix(eps)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Fe"].fix(eps)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sc"].fix(eps)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Y"].fix(225)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["La"].fix(5)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ce"].fix(50)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Pr"].fix(20)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Nd"].fix(75)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sm"].fix(50)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Gd"].fix(80)
m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Dy"].fix(60)

m.fs.solex_loading.mscontactor.aqueous_inlet_state[0].flow_vol.fix(50)

m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Al"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Ca"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Fe"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Sc"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Y"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["La"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Ce"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Pr"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Nd"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Sm"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Gd"].fix(eps)
m.fs.solex_loading.mscontactor.organic_inlet_state[0].conc_mass_comp["Dy"].fix(eps)

m.fs.solex_loading.mscontactor.organic_inlet_state[0].flow_vol.fix(50)

m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H2O"].fix(1e6)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["SO4"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["HSO4"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Al"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ca"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Fe"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sc"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Y"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["La"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ce"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Pr"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Nd"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sm"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Gd"].fix(eps)
m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Dy"].fix(eps)

m.fs.solex_stripping.mscontactor.aqueous_inlet_state[0].flow_vol.fix(50)

print(dof(m))

# Initialize and solve model

initializer = BlockTriangularizationInitializer(constraint_tolerance=1e-4)
initializer.initialize(m)

solver = SolverFactory("ipopt")
solver.options["bound_push"] = 1e-8
solver.options["mu_init"] = 1e-8
solver.solve(m, tee=True)

# Loading outlet

m.fs.solex_loading.mscontactor.organic[0,1].conc_mass_comp.pprint()

# Stripping outlet

m.fs.solex_stripping.mscontactor.aqueous[0,1].conc_mass_comp.pprint()