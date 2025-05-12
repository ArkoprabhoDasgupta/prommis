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
from prommis.solvent_extraction.membrane_channel_property_package import (
    MembraneSXChannelParameters,
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

m.fs.mem_prop = MembraneSXChannelParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = ReactionParameterTestBlock(property_package=m.fs.leach_soln)

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
    finite_elements=1,
    transformation_method="dae.finite_difference",
    transformation_scheme="BACKWARD",
    reaction_package=m.fs.reaxn,
)

m.fs.membrane_module.module_length.fix(1)
m.fs.membrane_module.tube_inner_radius.fix(0.1)
m.fs.membrane_module.tube_outer_radius.fix(0.2)
m.fs.membrane_module.shell_radius.fix(0.7)
m.fs.membrane_module.number_of_tubes.fix(5)

m.fs.membrane_module.feed_phase_inlet.flow_vol.fix(100)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp.fix(1e5)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.membrane_module.feed_phase_inlet.temperature.fix(303)
m.fs.membrane_module.feed_phase_inlet.pressure.fix(101325)

m.fs.membrane_module.strip_phase_inlet.flow_vol.fix(100)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp.fix(10)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.membrane_module.strip_phase_inlet.temperature.fix(303)
m.fs.membrane_module.strip_phase_inlet.pressure.fix(101325)

m.fs.membrane_module.feed_phase.properties[:, :].temperature.fix(303)
m.fs.membrane_module.feed_phase.properties[:, :].pressure.fix(101325)
m.fs.membrane_module.strip_phase.properties[:, :].temperature.fix(303)
m.fs.membrane_module.strip_phase.properties[:, :].pressure.fix(101325)

m.fs.membrane_module.feed_phase.mass_transfer_term.fix(-1)
m.fs.membrane_module.feed_phase.mass_transfer_term[:, :, "liquid", "H2O"].fix(0)
m.fs.membrane_module.strip_phase.mass_transfer_term.fix(1)
m.fs.membrane_module.strip_phase.mass_transfer_term[:, :, "liquid", "H2O"].fix(0)

# initializer = BlockTriangularizationInitializer()
# initializer.initialize(m)

print("Degrees of freedom:", dof(m))

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 10000
# solver.options["tol"] = 1e-5
results = solver.solve(m, tee=True)
