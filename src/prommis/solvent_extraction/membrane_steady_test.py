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
from prommis.solvent_extraction.membrane_solvent_extraction import (
    MembraneSolventExtraction,
)
from prommis.solvent_extraction.membrane_sx_transfer_package import (
    ReactionParameterTestBlock,
)
from prommis.leaching.leach_solution_properties import LeachSolutionParameters
import matplotlib.pyplot as plt

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.mem_prop = MembraneSXModuleParameters()
m.fs.leach_soln = LeachSolutionParameters()
# m.fs.reaxn = ReactionParameterTestBlock(property_package=m.fs.leach_soln)
m.fs.mem_prop.extractant_dosage = 3

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
    finite_elements=1,
    transformation_method="dae.finite_difference",
    transformation_scheme="BACKWARD",
    tube_inner_radius=0.1,
    tube_outer_radius=0.2,
    shell_radius=0.7,
    number_of_tubes=10,
    # reaction_package=m.fs.reaxn,
)

m.fs.membrane_module.module_length.fix(10.0)

m.fs.membrane_module.feed_phase_inlet.flow_vol.fix(100)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp.fix(10)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H"].fix(5.75)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "HSO4"].fix(1e3)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "SO4"].fix(100)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.membrane_module.feed_phase_inlet.temperature.fix(303)
m.fs.membrane_module.feed_phase_inlet.pressure.fix(101325)

m.fs.membrane_module.strip_phase_inlet.flow_vol.fix(100)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp.fix(1e-3)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "H"].fix(10.75)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "HSO4"].fix(1e3)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "SO4"].fix(100)
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


plt.plot(
    m.fs.membrane_module.feed_phase.length_domain,
    m.fs.membrane_module.feed_phase.properties[0, :].conc_mass_comp["Y"](),
)
plt.xlabel("Length, m")
plt.ylabel("Concentration, mg/L")
plt.title("Y feed concentration variation wrt length")
plt.figure()

plt.plot(
    m.fs.membrane_module.strip_phase.length_domain,
    m.fs.membrane_module.strip_phase.properties[0, :].conc_mass_comp["Y"](),
)
plt.xlabel("Length, m")
plt.ylabel("Concentration, mg/L")
plt.title("Y strip concentration variation wrt length")

percentage_recovery_feed = {}
REE_set = ["Y", "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy"]

for e in REE_set:
    percentage_recovery_feed[e] = [
        (
            1
            - (
                m.fs.membrane_module.feed_phase.properties[0, z].conc_mass_comp[e]()
                * m.fs.membrane_module.feed_phase.properties[0, z].flow_vol()
            )
            / (
                m.fs.membrane_module.feed_phase.properties[0, 0].conc_mass_comp[e]()
                * m.fs.membrane_module.feed_phase.properties[0, 0].flow_vol()
            )
        )
        * 100
        for z in m.fs.membrane_module.feed_phase.length_domain
    ]

percentage_recovery_strip = {}
REE_set = ["Y", "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy"]

for e in REE_set:
    percentage_recovery_strip[e] = [
        (
            (
                m.fs.membrane_module.strip_phase.properties[0, z].conc_mass_comp[e]()
                * m.fs.membrane_module.strip_phase.properties[0, z].flow_vol()
            )
            / (
                m.fs.membrane_module.strip_phase.properties[0, 0].conc_mass_comp[e]()
                * m.fs.membrane_module.strip_phase.properties[0, 0].flow_vol()
            )
            - 1
        )
        * 100
        for z in m.fs.membrane_module.strip_phase.length_domain
    ]
