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
from idaes.core.util import to_json, StoreSpec
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
m.fs.mem_prop.extractant_dosage = 10

m.fs.membrane_module = MembraneSolventExtraction(
    dynamic=False,
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
    finite_elements=4,
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

# m.discretizer = TransformationFactory("dae.finite_difference")
# m.discretizer.apply_to(m, nfe=len(time_set) - 1, wrt=m.fs.time, scheme="BACKWARD")

# Arcs connection

m.fs.feed_to_module = Arc(
    source=m.fs.feed_tank.outlet, destination=m.fs.membrane_module.feed_phase_inlet
)

m.fs.strip_to_module = Arc(
    source=m.fs.strip_tank.outlet, destination=m.fs.membrane_module.strip_phase_inlet
)

# m.fs.module_to_strip = Arc(
#     source=m.fs.membrane_module.strip_phase_outlet,
#     destination=m.fs.strip_tank.inlet,
# )

TransformationFactory("network.expand_arcs").apply_to(m)


# Define membrane module parameters

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

# Define feed tank inlet conditions

feed_flow_rate_values = {
    (0, 42): 15,
    (42, 72): 30,
    (72, 102): 45,
    (102, 112): 75,
    (112, 122): 150,
}

feed_concentration_values = {
    (0, 52): {
        "Al": 1338000,
        "Ca": 5720000,
        "Fe": 508344,
        "Y": 1111,
        "La": 1218,
        "Ce": 2506,
        "Pr": 331,
        "Nd": 1332,
        "Sm": 273,
        "Gd": 288,
        "Dy": 239,
        "Sc": 16.6,
    },
    (52, 82): {
        "Al": 1371000,
        "Ca": 5760000,
        "Fe": 468217,
        "Y": 1044,
        "La": 1185,
        "Ce": 2434,
        "Pr": 320,
        "Nd": 1309,
        "Sm": 269,
        "Gd": 274,
        "Dy": 226,
        "Sc": 16.8,
    },
    (82, 102): {
        "Al": 1275000,
        "Ca": 5562000,
        "Fe": 427389,
        "Y": 1017,
        "La": 1136,
        "Ce": 2323,
        "Pr": 306,
        "Nd": 1248,
        "Sm": 258,
        "Gd": 272,
        "Dy": 223,
        "Sc": 21.3,
    },
    (102, 122): {
        "Al": 1334000,
        "Ca": 5897000,
        "Fe": 448437,
        "Y": 1037,
        "La": 1145,
        "Ce": 2393,
        "Pr": 313,
        "Nd": 1268,
        "Sm": 261,
        "Gd": 274,
        "Dy": 221,
        "Sc": 18.9,
    },
}

for t in m.fs.time:
    if t == 0:
        m.fs.feed_tank.inlet.flow_vol[t].fix(
            feed_flow_rate_values[(0, 42)] * units.mL / units.min
        )
    else:
        for k, v in feed_flow_rate_values.items():
            if k[0] < t <= k[1]:
                m.fs.feed_tank.inlet.flow_vol[t].fix(v * units.mL / units.min)

    for k in feed_concentration_values.keys():
        for e, v in feed_concentration_values[k].items():
            if t == 0:
                m.fs.feed_tank.inlet.conc_mass_comp[t, e].fix(
                    feed_concentration_values[(0, 52)][e] * units.microgram / units.L
                )
            elif k[0] < t <= k[1]:
                m.fs.feed_tank.inlet.conc_mass_comp[t, e].fix(
                    v * units.microgram / units.L
                )

m.fs.feed_tank.inlet.conc_mass_comp[:, "H"].fix(10**-2 * 1 * units.gram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "H2O"].fix(1e6 * units.mg / units.L)
m.fs.feed_tank.control_volume.volume[:].fix(1 * units.L)

# Strip phase inlet conditions

strip_flow_rate_values = {
    (0, 12): 45,
    (12, 22): 75,
    (22, 42): 150,
    (42, 52): 45,
    (52, 63): 75,
    (63, 72): 150,
    (72, 82): 45,
    (82, 92): 75,
    (92, 122): 150,
}

for t in m.fs.time:
    if t == 0:
        m.fs.strip_tank.inlet.flow_vol[t].fix(45 * units.mL / units.min)
    else:
        for k, v in strip_flow_rate_values.items():
            if k[0] < t <= k[1]:
                m.fs.strip_tank.inlet.flow_vol[t].fix(v * units.mL / units.min)

m.fs.strip_tank.volume[:].fix(1 * units.L)

m.fs.strip_tank.inlet.conc_mass_comp[:, "Ca"].fix(578.1 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Al"].fix(1e-7 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Fe"].fix(92 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "La"].fix(0.211 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Ce"].fix(0.61 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Pr"].fix(0.029 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Nd"].fix(1e-7 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Sm"].fix(1e-7 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Gd"].fix(0.013 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Dy"].fix(0.046 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Y"].fix(0.102 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Sc"].fix(1e-7 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "H"].fix(1 * units.gram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "H2O"].fix(1e6 * units.mg / units.L)

print(dof(m))

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
solver.options["max_iter"] = 4000
results = solver.solve(scaled_model, tee=True)
# results = solver.solve(m, tee=True)

scaling.propagate_solution(scaled_model, m)

wts = StoreSpec.value()
to_json(m, fname="membrane_solvent_extraction_2.json", wts=wts, human_read=True)
