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

# assert 1 == 2

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

# Arcs connection

m.fs.feed_to_module = Arc(
    source=m.fs.feed_tank.outlet, destination=m.fs.membrane_module.feed_phase_inlet
)

m.fs.strip_to_module = Arc(
    source=m.fs.strip_tank.outlet, destination=m.fs.membrane_module.strip_phase_inlet
)


TransformationFactory("network.expand_arcs").apply_to(m)


# Define membrane module parameters

m.fs.membrane_module.module_length.fix(0.254 * 2 * units.m)
m.fs.membrane_module.eff[:, :, "Al"].fix(3.246e-4)
m.fs.membrane_module.eff[:, :, "Ca"].fix(1.765e-3)
m.fs.membrane_module.eff[:, :, "Ce"].fix(0.0305)
m.fs.membrane_module.eff[:, :, "Dy"].fix(1.156e-3)
m.fs.membrane_module.eff[:, :, "Fe"].fix(7.773e-4)
m.fs.membrane_module.eff[:, :, "Gd"].fix(0.013)
m.fs.membrane_module.eff[:, :, "La"].fix(0.032)
m.fs.membrane_module.eff[:, :, "Nd"].fix(0.059)
m.fs.membrane_module.eff[:, :, "Pr"].fix(0.105)
m.fs.membrane_module.eff[:, :, "Sm"].fix(0.15)
m.fs.membrane_module.eff[:, :, "Y"].fix(3e-6)
m.fs.membrane_module.eff[:, :, "Sc"].fix(1.5e-5)

# Define feed tank inlet conditions

feed_flow_rate_values = {
    (0, 10): 45,
    (10, 37): 15,
    (37, 64): 30,
    (64, 91): 45,
    (91, 97): 75,
    (97, 101): 150,
}

for t in m.fs.time:
    if t == 0:
        m.fs.feed_tank.inlet.flow_vol[t].fix(
            feed_flow_rate_values[(0, 10)] * units.mL / units.min
        )
    else:
        for k, v in feed_flow_rate_values.items():
            if k[0] < t <= k[1]:
                m.fs.feed_tank.inlet.flow_vol[t].fix(v * units.mL / units.min)

m.fs.feed_tank.inlet.conc_mass_comp[:, "Al"].fix(1351773 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Ca"].fix(6276510 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Fe"].fix(533032 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Sc"].fix(51.4 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "La"].fix(995 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Ce"].fix(1986 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Pr"].fix(275 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Nd"].fix(1104 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Sm"].fix(230 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Gd"].fix(231 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Dy"].fix(199 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Y"].fix(988 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "H"].fix(10**-2 * 1 * units.gram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "H2O"].fix(1e6 * units.mg / units.L)
m.fs.feed_tank.control_volume.volume[:].fix(4 * units.L)

# Strip phase inlet conditions

strip_flow_rate_values = {
    (0, 10): 150,
    (10, 19): 45,
    (19, 28): 75,
    (28, 37): 150,
    (37, 46): 45,
    (46, 55): 75,
    (55, 64): 150,
    (64, 73): 45,
    (73, 82): 75,
    (82, 101): 150,
}

for t in m.fs.time:
    if t == 0:
        m.fs.strip_tank.inlet.flow_vol[t].fix(150 * units.mL / units.min)
    else:
        for k, v in strip_flow_rate_values.items():
            if k[0] < t <= k[1]:
                m.fs.strip_tank.inlet.flow_vol[t].fix(v * units.mL / units.min)


m.fs.strip_tank.control_volume.volume[:].fix(0.72 * units.L)

# strip tank at t=0
m.fs.strip_tank.inlet.conc_mass_comp[:, "Ca"].fix(1e-8 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Al"].fix(860 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Fe"].fix(1e-8 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "La"].fix(0.12 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Ce"].fix(1.1 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Pr"].fix(1e-4 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Nd"].fix(0.72 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Sm"].fix(0.11 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Gd"].fix(0.1 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Dy"].fix(1e-8 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Y"].fix(1.4 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "Sc"].fix(2.59 * units.microgram / units.L)
m.fs.strip_tank.inlet.conc_mass_comp[:, "H"].fix(5 * units.gram / units.L)
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
    set_scaling_factor(
        m.fs.feed_tank.control_volume.material_balances[t, "H"],
        1e2,
    )
    for z in m.fs.membrane_module.feed_phase.length_domain:
        set_scaling_factor(
            m.fs.membrane_module.feed_phase.properties[t, z].conc_mol_comp["H"],
            1e2,
        )
        set_scaling_factor(
            m.fs.membrane_module.feed_phase.properties[t, z].pH_constraint,
            1e2,
        )
        set_scaling_factor(
            m.fs.membrane_module.feed_phase.properties[t, z].pH_phase,
            1e2,
        )
        set_scaling_factor(
            m.fs.membrane_module.feed_phase.material_flow_linking_constraints[
                t, z, "liquid", "H"
            ],
            1e2,
        )
        for e in m.fs.mem_prop.component_list:
            set_scaling_factor(
                m.fs.membrane_module.feed_phase.properties[t, z].conc_mol_comp[e],
                1e1,
            )
            set_scaling_factor(
                m.fs.membrane_module.feed_phase.properties[t, z].conc_mass_comp[e],
                1e2,
            )
            set_scaling_factor(
                m.fs.membrane_module.strip_phase.properties[t, z].conc_mol_comp[e],
                1e2,
            )

            for r in m.fs.membrane_module.r:
                set_scaling_factor(
                    m.fs.membrane_module.conc_mol_membrane_comp[t, z, r, e],
                    1e2,
                )


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

scaling.propagate_solution(scaled_model, m)

# scaling.propagate_solution(scaled_model, m)

to_json(m, fname="membrane_solvent_extraction_5M_HCl.json")
