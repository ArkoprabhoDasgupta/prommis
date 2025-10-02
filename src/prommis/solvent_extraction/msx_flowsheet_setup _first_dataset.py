from pyomo.environ import (
    ConcreteModel,
    units,
    Suffix,
    TransformationFactory,
    Var,
    minimize,
    log10,
)
from pyomo.network import Arc, SequentialDecomposition
from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    ControlVolume0DBlock,
)
from idaes.core.util import to_json
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.initialization import (
    SingleControlVolumeUnitInitializer,
    BlockTriangularizationInitializer,
)

from prommis.solvent_extraction.membrane_module_property_package import (
    MembraneSXModuleParameters,
)
from prommis.solvent_extraction.non_reactive_tank import NonReactiveTank
from prommis.solvent_extraction.membrane_channel_property_package import (
    MembraneSXChannelParameters,
)

import matplotlib.pyplot as plt
from prommis.solvent_extraction.membrane_solvent_extraction import (
    MembraneSolventExtraction,
    MembraneSolventExtractionInitializer,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from idaes.core.util import DiagnosticsToolbox

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.mem_prop = MembraneSXModuleParameters()
m.fs.mem_channel = MembraneSXChannelParameters()
# m.fs.leach_soln = LeachSolutionParameters()
m.fs.mem_prop.extractant_dosage = 10

m.fs.membrane_module = MembraneSolventExtraction(
    feed_phase={
        "property_package": m.fs.mem_channel,
        "flow_direction": FlowDirection.forward,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "energy_balance_type": EnergyBalanceType.none,
    },
    strip_phase={
        "property_package": m.fs.mem_channel,
        "flow_direction": FlowDirection.backward,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "energy_balance_type": EnergyBalanceType.none,
    },
    membrane_phase={
        "property_package": m.fs.mem_prop,
    },
    finite_elements=25,
    transformation_method="dae.finite_difference",
    transformation_scheme="BACKWARD",
    # collocation_points=2,
    tube_inner_radius=100 * units.micrometer,
    tube_outer_radius=150 * units.micrometer,
    shell_radius=(67 / 2) * units.mm,
    number_of_tubes=6300,
)

m.fs.feed_tank = NonReactiveTank(
    material_balance_type=MaterialBalanceType.componentTotal,
    dynamic=False,
    has_holdup=True,
    property_package=m.fs.mem_channel,
)

m.fs.strip_tank = NonReactiveTank(
    material_balance_type=MaterialBalanceType.componentTotal,
    dynamic=False,
    has_holdup=True,
    property_package=m.fs.mem_channel,
)


m.fs.membrane_module.module_length.fix(0.254 * 2 * units.m)
m.fs.membrane_module.eff["Al"].fix(3.246e-4)
m.fs.membrane_module.eff["Ca"].fix(1.765e-3)
m.fs.membrane_module.eff["Ce"].fix(0.0305)
m.fs.membrane_module.eff["Dy"].fix(1.156e-3)
m.fs.membrane_module.eff["Fe"].fix(7.773e-4)
m.fs.membrane_module.eff["Gd"].fix(0.013)
m.fs.membrane_module.eff["La"].fix(0.032)
m.fs.membrane_module.eff["Nd"].fix(0.059)
m.fs.membrane_module.eff["Pr"].fix(0.105)
m.fs.membrane_module.eff["Sm"].fix(0.15)
m.fs.membrane_module.eff["Y"].fix(3e-6)
m.fs.membrane_module.eff["Sc"].fix(1.5e-5)

print("check")

m.fs.feed_tank.inlet.flow_vol.fix(17 * units.mL / units.min)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Ca"].fix(1535976 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Al"].fix(310800 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Fe"].fix(131300 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "La"].fix(273 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Ce"].fix(573 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Pr"].fix(70 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Nd"].fix(315 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Sm"].fix(65.9 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Gd"].fix(67.5 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Dy"].fix(52.4 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Y"].fix(231 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "Sc"].fix(20 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "H"].fix(10**-2.06 * 1 * units.gram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[0, "H2O"].fix(1e6 * units.mg / units.L)
m.fs.feed_tank.volume.fix(10 * units.L)

m.fs.strip_tank.inlet.flow_vol.fix(54 * units.mL / units.min)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Ca"].fix(4081 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Al"].fix(1e-7 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Fe"].fix(515 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "La"].fix(1.89 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Ce"].fix(3.19 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Pr"].fix(0.398 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Nd"].fix(1.24 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Sm"].fix(0.15 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Gd"].fix(0.157 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Dy"].fix(0.119 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Y"].fix(0.6 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "Sc"].fix(0.649 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "H"].fix(3 * 1 * units.gram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[0, "H2O"].fix(1e6 * units.mg / units.L)
m.fs.strip_tank.volume.fix(10 * units.L)

m.fs.feed_to_module = Arc(
    source=m.fs.feed_tank.outlet, destination=m.fs.membrane_module.feed_phase_inlet
)

m.fs.strip_to_module = Arc(
    source=m.fs.strip_tank.outlet, destination=m.fs.membrane_module.strip_phase_inlet
)

TransformationFactory("network.expand_arcs").apply_to(m)

seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 5

# Using the SD tool
G = seq.create_graph(m)
heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
order = seq.calculation_order(G)

m.scaling_factor = Suffix(direction=Suffix.EXPORT)

for t in m.fs.time:
    for z in m.fs.membrane_module.feed_phase.length_domain:
        set_scaling_factor(
            m.fs.membrane_module.feed_phase.properties[t, z].conc_mol_comp["H"],
            1e4,
        )
        set_scaling_factor(
            m.fs.membrane_module.feed_phase.properties[t, z].pH_constraint,
            1e4,
        )
        # set_scaling_factor(
        #     m.fs.membrane_module.feed_phase.properties[t, z].pH_phase,
        #     1e2,
        # )
        set_scaling_factor(
            m.fs.membrane_module.feed_phase.material_flow_linking_constraints[
                t, z, "liquid", "H"
            ],
            1e4,
        )
        # for e in m.fs.mem_prop.component_list:
        #     set_scaling_factor(
        #         m.fs.membrane_module.feed_phase.properties[t, z].conc_mol_comp[e],
        #         1e2,
        #     )
        #     set_scaling_factor(
        #         m.fs.membrane_module.feed_phase.properties[t, z].conc_mass_comp[e],
        #         1e2,
        #     )
        #     set_scaling_factor(
        #         m.fs.membrane_module.strip_phase.properties[t, z].conc_mol_comp[e],
        #         1e2,
        #     )

        #     for r in m.fs.membrane_module.r:
        #         set_scaling_factor(
        #             m.fs.membrane_module.conc_mol_membrane_comp[t, z, r, e],
        #             1e2,
        #         )


scaling = TransformationFactory("core.scale_model")
scaled_model = scaling.create_using(m, rename=False)

tank_units = [scaled_model.fs.feed_tank, scaled_model.fs.strip_tank]
msx_units = [scaled_model.fs.membrane_module]


def function(unit):
    if unit in tank_units:
        BlockTriangularizationInitializer().initialize(unit)
    elif unit in msx_units:
        MembraneSolventExtractionInitializer().initialize(unit)


seq.run(scaled_model, function)

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 2000
results = solver.solve(scaled_model, tee=True)

scaling.propagate_solution(scaled_model, m)

to_json(m, fname="membrane_solvent_extraction.json")
