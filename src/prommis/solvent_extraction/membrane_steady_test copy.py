from pyomo.environ import (
    ConcreteModel,
    units,
)
import numpy as np
from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    SingleControlVolumeUnitInitializer,
)
from prommis.solvent_extraction.membrane_module_property_package import (
    MembraneSXModuleParameters,
)
from prommis.solvent_extraction.membrane_steady_state_model import (
    MembraneSolventExtraction,
)
from prommis.solvent_extraction.membrane_sx_transfer_package import (
    ReactionParameterTestBlock,
)
from prommis.leaching.leach_solution_properties import LeachSolutionParameters


m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.mem_prop = MembraneSXModuleParameters()
m.fs.leach_soln = LeachSolutionParameters()
# m.fs.reaxn = ReactionParameterTestBlock(property_package=m.fs.leach_soln)
m.fs.mem_prop.extractant_dosage = 1

m.fs.membrane_module = MembraneSolventExtraction(
    feed_phase={
        "property_package": m.fs.leach_soln,
        "flow_direction": FlowDirection.forward,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "energy_balance_type": EnergyBalanceType.none,
    },
    strip_phase={
        "property_package": m.fs.leach_soln,
        "flow_direction": FlowDirection.forward,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "energy_balance_type": EnergyBalanceType.none,
    },
    membrane_phase={
        "property_package": m.fs.mem_prop,
    },
    finite_elements=10,
    transformation_method="dae.finite_difference",
    transformation_scheme="BACKWARD",
    tube_inner_radius=100 * 1e-6,
    tube_outer_radius=150 * 1e-6,
    shell_radius=67e-3,
    number_of_tubes=6300,
    # reaction_package=m.fs.reaxn,
)

m.fs.membrane_module.module_length.fix(0.254)

m.fs.membrane_module.feed_phase_inlet.flow_vol.fix(6.3 * 60 / (6300 * 2820))
# m.fs.membrane_module.feed_phase_inlet.conc_mass_comp.fix(1e5)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Ca"].fix(664)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Al"].fix(3140)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Fe"].fix(109)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "La"].fix(1.73)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Ce"].fix(3.94)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Pr"].fix(0.54)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Nd"].fix(2.31)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Sm"].fix(0.566)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Gd"].fix(0.654)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Dy"].fix(0.64)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Y"].fix(3.65)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Sc"].fix(1)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Cl"].fix(1e-5)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H"].fix(0.186)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "HSO4"].fix(1e-5)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "SO4"].fix(8.938)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.membrane_module.feed_phase_inlet.temperature.fix(303)
m.fs.membrane_module.feed_phase_inlet.pressure.fix(101325)

m.fs.membrane_module.strip_phase_inlet.flow_vol.fix(1.5 * 60 / 2820)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp.fix(1e-5)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "H"].fix(0.5 * 1e3)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Cl"].fix(0.5 * 35.5e3)
# m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "HSO4"].fix(1e3)
# m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "SO4"].fix(100)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.membrane_module.strip_phase_inlet.temperature.fix(303)
m.fs.membrane_module.strip_phase_inlet.pressure.fix(101325)

m.fs.membrane_module.feed_phase.properties[:, :].temperature.fix(303)
m.fs.membrane_module.feed_phase.properties[:, :].pressure.fix(101325)
m.fs.membrane_module.strip_phase.properties[:, :].temperature.fix(303)
m.fs.membrane_module.strip_phase.properties[:, :].pressure.fix(101325)

# m.fs.membrane_module.feed_phase.mass_transfer_term.fix(-1)
# m.fs.membrane_module.feed_phase.mass_transfer_term[:, :, "liquid", "H2O"].fix(0)
# m.fs.membrane_module.strip_phase.mass_transfer_term.fix(1)
# m.fs.membrane_module.strip_phase.mass_transfer_term[:, :, "liquid", "H2O"].fix(0)

# initializer = BlockTriangularizationInitializer()
# initializer.initialize(m)

print("Degrees of freedom:", dof(m))

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 10000
# # solver.options["tol"] = 1e-5
results = solver.solve(m, tee=True)
