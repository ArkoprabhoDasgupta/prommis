#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import ConcreteModel, units, TransformationFactory, RangeSet

from idaes.core import FlowDirection, FlowsheetBlock, MaterialBalanceType
from idaes.core.util import to_json

from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.solvers import get_solver
from pyomo.network import Arc, SequentialDecomposition
from prommis.solvent_extraction.leach_solution_properties_optimization import (
    LeachSolutionParameters,
)
from prommis.solvent_extraction.ree_og_distribution_new_optimization import (
    REESolExOgParameters,
)
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
    MixerSettlerExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction import (
    SolventExtraction,
    SolventExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package_new_modified_for_flowsheet_optimization import (
    SolventExtractionReactions,
)
from prommis.solvent_extraction.neutralization_tank import NeutralizationTank
from idaes.models.unit_models.mixer import (
    Mixer,
    MixerInitializer,
    MixingType,
    MomentumMixingType,
)
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
)


m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()


# define stages
dosage = 8

# define scrub sx
m.fs.scrub_sx = SolventExtraction(
    number_of_finite_elements=3,
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
)


# define inlet conditions

# define scrub sx aqueous inlet
m.fs.scrub_sx.aqueous_inlet.flow_vol.fix(60.01)
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "Cl"]:
        m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, e].fix(1e-9)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H"].fix(6 * units.g / units.L)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(6 * 35.5 * units.g / units.L)

m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Ce_o"].fix(0.635)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "DEHPA"].fix(77952.034)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Dy_o"].fix(0.083)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(5.573)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820000)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "La_o"].fix(0.2277)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Nd_o"].fix(0.2705)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Pr_o"].fix(0.1607)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Sm_o"].fix(0.059)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Y_o"].fix(0.319)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Gd_o"].fix(0.3559)
m.fs.scrub_sx.organic_inlet.flow_vol.fix(65.01)

m.fs.scrub_sx.mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.scrub_sx.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
m.fs.scrub_sx.mscontactor.organic[:, :].temperature.fix(305.15 * units.K)


solver = get_solver("ipopt_v2")
solver.options["halt_on_ampl_error"] = "yes"
solver.options["max_iter"] = 1000

print(dof(m))


sx_initializer = SolventExtractionInitializer()
scrub_sx_units = m.fs.scrub_sx

sx_initializer.initialize(scrub_sx_units)

results = solver.solve(m, tee=True)
