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
load_number_of_stages = 4
load_stage_list = RangeSet(1, load_number_of_stages)
load_interstage_list = RangeSet(1, load_number_of_stages - 1)

strip_number_of_stages = 3
strip_stage_list = RangeSet(1, strip_number_of_stages)
strip_interstage_list = RangeSet(1, strip_number_of_stages - 1)

m.fs.reaxn.extractant_dosage = dosage

# define load sx
m.fs.load_sx = MixerSettlerExtraction(
    load_stage_list,
    number_of_stages=1,
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
    settler_finite_elements=2,
)

# define scrub sx
m.fs.scrub_sx = MixerSettlerExtraction(
    number_of_stages=1,
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
    settler_finite_elements=2,
)

# define strip sx
m.fs.strip_sx = MixerSettlerExtraction(
    strip_stage_list,
    number_of_stages=1,
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
    settler_finite_elements=2,
)


# define neutralization tank
m.fs.aq_feed_neutral = NeutralizationTank(property_package=m.fs.leach_soln)

# define organic interstage mixer
m.fs.org_inter_mixer = Mixer(
    load_interstage_list,
    property_package=m.fs.prop_o,
    num_inlets=2,
    inlet_list=["sx", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

# define aqueous interstage mixer
m.fs.aq_inter_mixer = Mixer(
    strip_interstage_list,
    property_package=m.fs.leach_soln,
    num_inlets=2,
    inlet_list=["sx", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

# define arcs

# neutral tank to load sx
m.fs.neutral_to_sx = Arc(
    source=m.fs.aq_feed_neutral.outlet, destination=m.fs.load_sx[1].aqueous_inlet
)

for i in load_stage_list:

    if i != 1:
        # loading sx aqueous phases
        m.add_component(
            f"load_aqueous_sx_{i-1}_to_{i}",
            Arc(
                source=m.fs.load_sx[i - 1].aqueous_outlet,
                destination=m.fs.load_sx[i].aqueous_inlet,
            ),
        )
        # loading sx organic to interstage
        m.add_component(
            f"load_organic_sx_{i}_to_interstage_{i-1}",
            Arc(
                source=m.fs.load_sx[i].organic_outlet,
                destination=m.fs.org_inter_mixer[i - 1].sx,
            ),
        )
    # interstage to loading organic
    if i != load_number_of_stages:
        m.add_component(
            f"load_interstage_organic_{i}_to_sx_organic_{i}",
            Arc(
                source=m.fs.org_inter_mixer[i].outlet,
                destination=m.fs.load_sx[i].organic_inlet,
            ),
        )

for i in strip_stage_list:

    if i != strip_number_of_stages:

        # stripping sx mixer to aqueous sx
        m.add_component(
            f"strip_aqueous_sx_{i}_to_mixer_{i}",
            Arc(
                source=m.fs.strip_sx[i].aqueous_outlet,
                destination=m.fs.aq_inter_mixer[i].sx,
            ),
        )
        # strip mixer aq to sx aq inlet
        m.add_component(
            f"strip_aqueous_mixer_{i}_to_sx_{i+1}",
            Arc(
                source=m.fs.aq_inter_mixer[i].outlet,
                destination=m.fs.strip_sx[i + 1].aqueous_inlet,
            ),
        )
        # stripping organic
        m.add_component(
            f"strip_organic_sx_{i+1}_to_{i}",
            Arc(
                source=m.fs.strip_sx[i + 1].organic_outlet,
                destination=m.fs.strip_sx[i].organic_inlet,
            ),
        )

m.organic_load_to_scrub = Arc(
    source=m.fs.load_sx[1].organic_outlet,
    destination=m.fs.scrub_sx.organic_inlet,
)

m.organic_scrub_to_strip = Arc(
    source=m.fs.scrub_sx.organic_outlet,
    destination=m.fs.strip_sx[strip_number_of_stages].organic_inlet,
)


TransformationFactory("network.expand_arcs").apply_to(m)


# define inlet conditions

# define neutral tank conditions
pH = 1.3
m.fs.aq_feed_neutral.inlet.flow_vol.fix(60.01)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ascorbic"].fix(1.44 * units.g / units.L)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Fe"].fix(138.27)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "La"].fix(2.09)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ce"].fix(5)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Pr"].fix(0.73)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Nd"].fix(2.10)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Sm"].fix(0.236)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Gd"].fix(0.56)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Dy"].fix(0.09)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Y"].fix(0.346)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "H"].fix(10 ** (-pH) * units.g / units.L)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "SO4"].fix(
    10 ** (-pH) * 48 * units.g / units.L
)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "HSO4"].fix(1e-4)
m.fs.aq_feed_neutral.base_flowrate[0].fix(1)
m.fs.aq_feed_neutral.base_concentration[0].fix(0.1)

m.fs.aq_feed_neutral.control_volume.properties_out[0.0].pressure.fix(101235)
m.fs.aq_feed_neutral.control_volume.properties_out[0.0].temperature.fix(303.5)


# define load sx organic inlet
m.fs.load_sx[load_number_of_stages].organic_inlet.conc_mass_comp[0, "Kerosene"].fix(
    820e3
)
m.fs.load_sx[load_number_of_stages].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
    975.8e3 * dosage / 100
)
for e in m.fs.prop_o.component_list:
    if e not in ["Kerosene", "DEHPA"]:
        m.fs.load_sx[load_number_of_stages].organic_inlet.conc_mass_comp[0, e].fix(1e-7)
m.fs.load_sx[load_number_of_stages].organic_inlet.flow_vol.fix(62.01)

m.fs.load_sx[:].mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.load_sx[:].mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.load_sx[:].mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.load_sx[:].organic_settler[:].unit.area.fix(1)
m.fs.load_sx[:].aqueous_settler[:].unit.area.fix(1)
m.fs.load_sx[:].aqueous_settler[:].unit.length.fix(1)
m.fs.load_sx[:].organic_settler[:].unit.length.fix(1)


# define organic interstage inlet
for i in load_interstage_list:
    for e in m.fs.prop_o.component_list:
        if e not in ["Kerosene", "DEHPA"]:
            m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-7)
    m.fs.org_inter_mixer[i].feed.flow_vol.fix(1)
    m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "Kerosene"].fix(820e3)
    m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "DEHPA"].fix(
        975.8e3 * dosage / 100
    )  # dv

m.fs.org_inter_mixer[:].mixed_state[0.0].pressure.fix(101235)
m.fs.org_inter_mixer[:].mixed_state[0.0].temperature.fix(303.5)

# define aqueous interstage inlet
for i in strip_interstage_list:
    for e in m.fs.leach_soln.component_list:
        if e not in ["H2O", "H", "Cl"]:
            m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-9)
    m.fs.aq_inter_mixer[i].feed.flow_vol.fix(0.5)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H2O"].fix(1e6)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].fix(1 * units.g / units.L)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "Cl"].fix(
        1 * 35.5 * units.g / units.L
    )

m.fs.aq_inter_mixer[:].mixed_state[0.0].pressure.fix(101235)
m.fs.aq_inter_mixer[:].mixed_state[0.0].temperature.fix(303.5)

# define scrub sx aqueous inlet
m.fs.scrub_sx.aqueous_inlet.flow_vol.fix(60.01)
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "Cl"]:
        m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, e].fix(1e-9)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H"].fix(0.1 * units.g / units.L)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(0.1 * 35.5 * units.g / units.L)

m.fs.scrub_sx.mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.scrub_sx.organic_settler[:].unit.area.fix(1)
m.fs.scrub_sx.aqueous_settler[:].unit.area.fix(1)
m.fs.scrub_sx.aqueous_settler[:].unit.length.fix(1)
m.fs.scrub_sx.organic_settler[:].unit.length.fix(1)
m.fs.scrub_sx.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
m.fs.scrub_sx.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(305.15 * units.K)


# define strip sx aqueous inlet
m.fs.strip_sx[1].aqueous_inlet.flow_vol.fix(60.01)
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "Cl"]:
        m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, e].fix(1e-9)
m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "H"].fix(1 * units.g / units.L)
m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1 * 35.5 * units.g / units.L)
m.fs.strip_sx[:].mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.strip_sx[:].organic_settler[:].unit.area.fix(1)
m.fs.strip_sx[:].aqueous_settler[:].unit.area.fix(1)
m.fs.strip_sx[:].aqueous_settler[:].unit.length.fix(1)
m.fs.strip_sx[:].organic_settler[:].unit.length.fix(1)
m.fs.strip_sx[:].mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.strip_sx[:].mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
    305.15 * units.K
)

sx_initializer = MixerSettlerExtractionInitializer()
load_sx_units = [m.fs.load_sx[s] for s in load_stage_list]
scrub_sx_units = m.fs.scrub_sx
strip_sx_units = [m.fs.strip_sx[s] for s in strip_stage_list]

mixer_initializer = MixerInitializer()
load_interstage_mixer = [m.fs.org_inter_mixer[s] for s in load_interstage_list]
strip_interstage_mixer = [m.fs.aq_inter_mixer[s] for s in strip_interstage_list]

tank_initializer = BlockTriangularizationInitializer()
neutral_tank = m.fs.aq_feed_neutral


solver = get_solver("ipopt_v2")
solver.options["halt_on_ampl_error"] = "yes"
solver.options["max_iter"] = 1000

# identify tear streams by sequential decomposition
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 3

G = seq.create_graph(m)
tear_arcs = seq.tear_set_arcs(G, solver=solver)

print("Tear streams identified by sequential decomposition:")
for arc in tear_arcs:
    print(arc.name)


# def function(unit):
#     if unit in load_sx_units:
#         print(dof(unit))
#         print(f"Initializing {unit}")
#         sx_initializer.initialize(unit)
#     elif unit in scrub_sx_units:
#         print(f"Initializing {unit}")
#         sx_initializer.initialize(unit)
#     elif unit in strip_sx_units:
#         print(f"Initializing {unit}")
#         sx_initializer.initialize(unit)
#     elif unit in load_interstage_mixer:
#         print(dof(unit))
#         print(f"Initializing {unit}")
#         mixer_initializer.initialize(unit)
#     elif unit in strip_interstage_mixer:
#         print(dof(unit))
#         print(f"Initializing {unit}")
#         mixer_initializer.initialize(unit)
#     else:
#         print(f"Initializing {unit}")
#         tank_initializer.initialize(unit)


# seq.run(m, function)


print(dof(m))

# # initializer = m.fs.load_sx.default_initializer()
# # initializer.initialize(m.fs.load_sx)


results = solver.solve(m, tee=True)

# percentage_recovery = {}
# for e in m.fs.leach_soln.component_list:
#     if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic"]:
#         percentage_recovery[e] = [
#             (
#                 (
#                     m.fs.load_sx[1].organic_outlet.conc_mass_comp[
#                         t, f"{e}_o"
#                     ]()
#                     * m.fs.load_sx[1].organic_outlet.flow_vol[t]()
#                     - m.fs.load_sx[
#                         load_number_of_stages
#                     ].organic_inlet.conc_mass_comp[t, f"{e}_o"]()
#                     * m.fs.load_sx[load_number_of_stages].organic_inlet.flow_vol[
#                         t
#                     ]()
#                 )
#                 / (
#                     m.fs.aq_feed_neutral.inlet.conc_mass_comp[t, e]()
#                     * m.fs.aq_feed_neutral.inlet.flow_vol[t]()
#                 )
#             )
#             * 100
#             for t in m.fs.time
#         ]


# percentage_recovery["tree"] = [
#     (
#         (
#             sum(
#                 m.fs.load_sx[1].organic_outlet.conc_mass_comp[t, f"{e}_o"]()
#                 * m.fs.load_sx[1].organic_outlet.flow_vol[t]()
#                 for e in m.fs.leach_soln.component_list
#                 if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
#             )
#             - sum(
#                 m.fs.load_sx[load_number_of_stages].organic_inlet.conc_mass_comp[
#                     t, f"{e}_o"
#                 ]()
#                 * m.fs.load_sx[load_number_of_stages].organic_inlet.flow_vol[t]()
#                 for e in m.fs.leach_soln.component_list
#                 if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
#             )
#         )
#         / sum(
#             m.fs.aq_feed_neutral.inlet.conc_mass_comp[t, e]()
#             * m.fs.aq_feed_neutral.inlet.flow_vol[t]()
#             for e in m.fs.leach_soln.component_list
#             if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
#         )
#     )
#     * 100
#     for t in m.fs.time
# ]

# print(percentage_recovery)
