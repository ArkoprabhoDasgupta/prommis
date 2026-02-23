from matplotlib.category import _log
from pyomo.environ import (
    ConcreteModel,
    units,
    RangeSet,
    TransformationFactory,
    Var,
    maximize,
    Param,
    value,
    minimize,
    Constraint,
    Suffix,
)
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
)
from idaes.core.util import to_json
from idaes.core.util.model_diagnostics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.initialization import (
    SingleControlVolumeUnitInitializer,
    BlockTriangularizationInitializer,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution_new import REESolExOgParameters
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
    MixerSettlerExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package_new_modified import (
    SolventExtractionReactions,
)
from prommis.solvent_extraction.neutralization_tank import NeutralizationTank
from idaes.models.unit_models.mixer import (
    Mixer,
    MixerInitializer,
    MixingType,
    MomentumMixingType,
)

# Create the model

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()

# Define the number of stages in each section

loading_stages = 4
load_stage_list = RangeSet(1, loading_stages)
load_interstage_list = RangeSet(1, loading_stages - 1)

strip_stages = 3
strip_stage_list = RangeSet(1, strip_stages)

scrub_stages = 1
scrub_stage_list = RangeSet(1, scrub_stages)

# Define the sx units

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
    settler_finite_elements=4,
)

# m.fs.scrub_sx = MixerSettlerExtraction(
#     number_of_stages=1,
#     aqueous_stream={
#         "property_package": m.fs.leach_soln,
#         "flow_direction": FlowDirection.forward,
#         "has_energy_balance": False,
#         "has_pressure_balance": False,
#     },
#     organic_stream={
#         "property_package": m.fs.prop_o,
#         "flow_direction": FlowDirection.backward,
#         "has_energy_balance": False,
#         "has_pressure_balance": False,
#     },
#     heterogeneous_reaction_package=m.fs.reaxn,
#     has_holdup=True,
#     settler_transformation_method="dae.finite_difference",
#     settler_transformation_scheme="BACKWARD",
#     settler_finite_elements=4,
# )

# m.fs.strip_sx = MixerSettlerExtraction(
#     strip_stage_list,
#     number_of_stages=1,
#     aqueous_stream={
#         "property_package": m.fs.leach_soln,
#         "flow_direction": FlowDirection.forward,
#         "has_energy_balance": False,
#         "has_pressure_balance": False,
#     },
#     organic_stream={
#         "property_package": m.fs.prop_o,
#         "flow_direction": FlowDirection.backward,
#         "has_energy_balance": False,
#         "has_pressure_balance": False,
#     },
#     heterogeneous_reaction_package=m.fs.reaxn,
#     has_holdup=True,
#     settler_transformation_method="dae.finite_difference",
#     settler_transformation_scheme="BACKWARD",
#     settler_finite_elements=4,
# )

# define interstage mixers

m.fs.org_inter_mixer = Mixer(
    load_interstage_list,
    property_package=m.fs.prop_o,
    num_inlets=2,
    inlet_list=["sx", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

# m.fs.aq_inter_mixer = Mixer(
#     strip_stage_list,
#     property_package=m.fs.leach_soln,
#     num_inlets=2,
#     inlet_list=["sx", "feed"],
#     material_balance_type=MaterialBalanceType.componentTotal,
#     energy_mixing_type=MixingType.none,
#     momentum_mixing_type=MomentumMixingType.none,
# )

# add neutral tanks

m.fs.aq_feed_neutral = NeutralizationTank(property_package=m.fs.leach_soln)

# add arcs

# loading arcs

# neutral tank to sx loading aqueous inlet
m.neutral_to_load_aq_sx = Arc(
    source=m.fs.aq_feed_neutral.outlet, destination=m.fs.load_sx[1].aqueous_inlet
)

for i in load_stage_list:
    if i != 1:
        # loading sx aqueous phases
        m.add_component(
            f"load_aq_sx_{i-1}_to_{i}",
            Arc(
                source=m.fs.load_sx[i - 1].aqueous_outlet,
                destination=m.fs.load_sx[i].aqueous_inlet,
            ),
        )
        # loading sx organic to interstage
        m.add_component(
            f"load_org_sx_{i}_to_load_org_inter_{i-1}",
            Arc(
                source=m.fs.load_sx[i].organic_outlet,
                destination=m.fs.org_inter_mixer[i - 1].sx,
            ),
        )
    # interstage to loading organic
    if i != loading_stages:
        m.add_component(
            f"load_org_inter_{i}_to_load_org_sx_{i}",
            Arc(
                source=m.fs.org_inter_mixer[i].outlet,
                destination=m.fs.load_sx[i].organic_inlet,
            ),
        )

# # stripping arcs

# for i in strip_stage_list:
#     # stripping aqueous mixer to aqueous sx
#     m.add_component(
#         f"strip_aq_inter_{i}_to_strip_aq_sx_{i}",
#         Arc(
#             source=m.fs.aq_inter_mixer[i].outlet,
#             destination=m.fs.strip_sx[i].aqueous_inlet,
#         ),
#     )
#     if i != strip_stages:
#         # stripping sx mixer to aqueous sx
#         m.add_component(
#             f"strip_aq_sx_{i}_to_strip_aq_inter_{i+1}",
#             Arc(
#                 source=m.fs.strip_sx[i].aqueous_outlet,
#                 destination=m.fs.aq_inter_mixer[i + 1].sx,
#             ),
#         )
#         # stripping organic
#         m.add_component(
#             f"strip_org_sx_{i+1}_to_strip_organic_sx_{i}",
#             Arc(
#                 source=m.fs.strip_sx[i + 1].organic_outlet,
#                 destination=m.fs.strip_sx[i].organic_inlet,
#             ),
#         )

# # connect organic loading and scrubbing and stripping

# m.org_load_out_to_scrub_in = Arc(
#     source=m.fs.load_sx[1].organic_outlet,
#     destination=m.fs.scrub_sx.organic_inlet,
# )

# m.org_scrub_out_to_strip_in = Arc(
#     source=m.fs.scrub_sx.organic_outlet,
#     destination=m.fs.strip_sx[strip_stages].organic_inlet,
# )

TransformationFactory("network.expand_arcs").apply_to(m)

# set inlets

# neutralization tank input

pH_load = 0.6
m.fs.aq_feed_neutral.inlet.flow_vol.fix(62.01)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Al"].fix(137.27)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ca"].fix(25.78)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Fe"].fix(138.27)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Sc"].fix(0.277)
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
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "H"].fix(
    10 ** (-pH_load) * units.gram / units.L
)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "SO4"].fix(
    10 ** (-pH_load) * 48 * units.gram / units.L
)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "HSO4"].fix(1e-4)
m.fs.aq_feed_neutral.base_flowrate[0].fix(1)
m.fs.aq_feed_neutral.base_concentration[0].fix(0.1)  # dv

# loading sx organic inlet

for e in m.fs.prop_o.component_list:
    if e not in ["Kerosene", "DEHPA"]:
        m.fs.load_sx[loading_stages].organic_inlet.conc_mass_comp[0, e].fix(1e-9)
m.fs.load_sx[loading_stages].organic_inlet.flow_vol.fix(62.01)
m.fs.load_sx[loading_stages].organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)
dosage = 8  # in vol %
m.fs.load_sx[loading_stages].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
    975.8e3 * dosage / 100
)

# loading organic interstage addition

for i in load_interstage_list:
    for e in m.fs.prop_o.component_list:
        if e not in ["Kerosene", "DEHPA"]:
            m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-9)
    m.fs.org_inter_mixer[i].feed.flow_vol.fix(2)
    m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "Kerosene"].fix(820e3)
    dosage = 10  # in vol %
    m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "DEHPA"].fix(
        975.8e3 * dosage / 100
    )  # dv

# # stripping aqueous interstage addition

# pH_strip = 0.5
# for i in strip_stage_list:
#     for e in m.fs.leach_soln.component_list:
#         if e not in ["H2O", "HSO4", "SO4", "H"]:
#             m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-9)
#     m.fs.aq_inter_mixer[i].feed.flow_vol.fix(2)
#     m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H2O"].fix(1e6)
#     m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].fix(
#         10 ** (-pH_strip) * 2 * units.gram / units.L
#     )  # dv
#     m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "SO4"].fix(1e-7)
#     m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "HSO4"].fix(1e-7)

# # aqueous strip sx inlet

# pH_strip = 0.4

# for e in m.fs.leach_soln.component_list:
#     if e not in ["H2O", "HSO4", "SO4", "H"]:
#         m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, e].fix(1e-7)
# m.fs.aq_inter_mixer[1].sx.flow_vol.fix(62.01)
# m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "H2O"].fix(1e6)
# m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "H"].fix(
#     1 * units.gram / units.L
# )  # maybe dv
# m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "SO4"].fix(1e-7)
# m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "HSO4"].fix(1e-7)

# # aqueous scrub sx feed

# for e in m.fs.leach_soln.component_list:
#     if e not in ["H2O", "HSO4", "SO4", "H"]:
#         m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, e].fix(1e-9)
# m.fs.scrub_sx.aqueous_inlet.flow_vol.fix(1)
# m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
# m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H"].fix(10**-0.3 * units.gram / units.L)
# m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(1e-7)
# m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-7)


# fix parameters

# fix mixer settler volumes, areas, lengths, and temperatures

m.fs.load_sx[:].mixer[:].unit.mscontactor.volume[:].fix(1e-2 * units.m**3)
m.fs.load_sx[:].organic_settler[:].unit.area.fix(1e-2)
m.fs.load_sx[:].aqueous_settler[:].unit.area.fix(1e-2)
m.fs.load_sx[:].aqueous_settler[:].unit.length.fix(1e-2)
m.fs.load_sx[:].organic_settler[:].unit.length.fix(1e-2)
m.fs.load_sx[:].mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.load_sx[:].mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
    305.15 * units.K
)

# m.fs.strip_sx[:].mixer[:].unit.mscontactor.volume[:].fix(1e-2 * units.m**3)
# m.fs.strip_sx[:].organic_settler[:].unit.area.fix(1e-2)
# m.fs.strip_sx[:].aqueous_settler[:].unit.area.fix(1e-2)
# m.fs.strip_sx[:].aqueous_settler[:].unit.length.fix(1e-2)
# m.fs.strip_sx[:].organic_settler[:].unit.length.fix(1e-2)
# m.fs.strip_sx[:].mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
#     305.15 * units.K
# )
# m.fs.strip_sx[:].mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
#     305.15 * units.K
# )

# m.fs.scrub_sx.mixer[:].unit.mscontactor.volume[:].fix(1e-2 * units.m**3)
# m.fs.scrub_sx.organic_settler[:].unit.area.fix(1e-2)
# m.fs.scrub_sx.aqueous_settler[:].unit.area.fix(1e-2)
# m.fs.scrub_sx.aqueous_settler[:].unit.length.fix(1e-2)
# m.fs.scrub_sx.organic_settler[:].unit.length.fix(1e-2)
# m.fs.scrub_sx.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
# m.fs.scrub_sx.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(305.15 * units.K)

# fix neutral tank volumes and temperatures

m.fs.aq_feed_neutral.control_volume.properties_out[0.0].temperature.fix(305)
m.fs.aq_feed_neutral.control_volume.properties_out[0.0].pressure.fix(101325)

# fix interstage mixer

m.fs.org_inter_mixer[:].mixed_state[0.0].temperature.fix(305)
# m.fs.aq_inter_mixer[:].mixed_state[0.0].temperature.fix(305)
m.fs.org_inter_mixer[:].mixed_state[0.0].pressure.fix(101325)
# m.fs.aq_inter_mixer[:].mixed_state[0.0].pressure.fix(101325)

m.fs.load_sx[1].mixer[1].unit.mscontactor.heterogeneous_reactions[
    0.0, 1
].ascorbic_dosage.fix(0)
# m.fs.strip_sx[:].mixer[1].unit.mscontactor.heterogeneous_reactions[
#     0.0, 1
# ].ascorbic_dosage.fix(0)
# m.fs.scrub_sx.mixer[1].unit.mscontactor.heterogeneous_reactions[
#     0.0, 1
# ].ascorbic_dosage.fix(0)


@m.Constraint(load_stage_list)
def ascorbic_acid_constraint(m, s):
    if s == 1:
        return Constraint.Skip
    else:
        return (
            m.fs.load_sx[s]
            .mixer[1]
            .unit.mscontactor.heterogeneous_reactions[0.0, 1]
            .ascorbic_dosage
            == m.fs.load_sx[s - 1]
            .mixer[1]
            .unit.mscontactor.heterogeneous_reactions[0.0, 1]
            .ascorbic_dosage
        )


# assert 1==2

# m.scaling_factor = Suffix(direction=Suffix.EXPORT)

# # set_scaling_factor(m.fs.aq_inter_mixer[1].feed_state[0.0].pH_constraint['liquid'], 1e3)
# for s in strip_stage_list:
#     if s != 2:
#         set_scaling_factor(
#             m.fs.aq_inter_mixer[s].feed_state[0.0].pH_constraint["liquid"], 1e3
#         )
#     else:
#         set_scaling_factor(
#             m.fs.aq_inter_mixer[s].sx_state[0.0].pH_constraint["liquid"], 1e3
#         )
#     for e in REE_list:
#         set_scaling_factor(
#             m.fs.strip_sx[s].mixer[1].unit.distribution_extent_constraint[0, 1, e], 1
#         )
#     set_scaling_factor(
#         m.fs.strip_sx[s].mixer[1].unit.distribution_extent_constraint[0, 1, "Fe"], 1
#     )
#     set_scaling_factor(
#         m.fs.strip_sx[s].mixer[1].unit.distribution_extent_constraint[0, 1, "Al"], 1
#     )
# for s in load_stage_list:
#     for e in REE_list:
#         set_scaling_factor(
#             m.fs.load_sx[s].mixer[1].unit.distribution_extent_constraint[0, 1, e], 1
#         )


# scaling = TransformationFactory("core.scale_model")
# scaled_model = scaling.create_using(m, rename=False)

print(degrees_of_freedom(m))

# assert 1 == 2
# seq = SequentialDecomposition()
# seq.options.select_tear_method = "heuristic"
# seq.options.tear_method = "Wegstein"
# seq.options.iterLim = 3

# # Using the SD tool
# G = seq.create_graph(m)
# heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
# order = seq.calculation_order(G)

# for o in heuristic_tear_set:
#     print(o.name)

# Assuming 'model' is your Pyomo model with Network components
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
graph = seq.create_graph(m)

# Identify tear sets automatically
tears = seq.tear_set_arcs(graph, solver="ipopt_v2")

# Print the tear streams
print("Tear Streams Identified:")
for edge in tears:
    print(edge)  # This will show the source/destination port information


# seq = SequentialDecomposition()
# seq.options.tear_method = "Direct"
# seq.options.tear_solver = "ipopt_v2"
# seq.options.iterLim = 1
# # seq.options.tear_set = [
# #     m.load_aq_sx_1_to_2,
# #     m.load_aq_sx_2_to_3,
# #     m.load_aq_sx_3_to_4,
# #     m.strip_aq_sx_1_to_strip_aq_inter_2,
# #     m.strip_aq_sx_2_to_strip_aq_inter_3,
# # ]

# assert 1 == 2
tear_guesses1 = {
    "flow_vol": {0: 62.01},
    "conc_mass_comp": {
        (0, "Al"): 157.27,
        (0, "Ca"): 25.78,
        (0, "Ce"): 5,
        (0, "Cl"): 1e-7,
        (0, "Dy"): 0.09,
        (0, "Fe"): 138.27,
        (0, "Gd"): 0.56,
        (0, "H"): 10 ** (3 - 0.6),
        (0, "H2O"): 1000000,
        (0, "HSO4"): 1e-4,
        (0, "La"): 2.09,
        (0, "Nd"): 2.1,
        (0, "Pr"): 0.73,
        (0, "SO4"): 10 ** (3 - 0.6) * 48,
        (0, "Sc"): 0.277,
        (0, "Sm"): 0.236,
        (0, "Y"): 0.346,
    },
}

# tear_guesses2 = {
#     "flow_vol": {0: 62.01},
#     "conc_mass_comp": {
#         (0, "Al"): 1e-4,
#         (0, "Ca"): 1e-4,
#         (0, "Ce"): 1e-4,
#         (0, "Cl"): 1e-4,
#         (0, "Dy"): 1e-4,
#         (0, "Fe"): 1e-4,
#         (0, "Gd"): 0.56,
#         (0, "H"): 1000,
#         (0, "H2O"): 1000000,
#         (0, "HSO4"): 1e-4,
#         (0, "La"): 1e-4,
#         (0, "Nd"): 1e-4,
#         (0, "Pr"): 1e-4,
#         (0, "SO4"): 1e-4,
#         (0, "Sc"): 1e-4,
#         (0, "Sm"): 1e-4,
#         (0, "Y"): 1e-4,
#     },
# }

# seq.set_guesses_for(m.fs.load_sx[2].aqueous_inlet, tear_guesses1)
# seq.set_guesses_for(m.fs.load_sx[3].aqueous_inlet, tear_guesses1)
# seq.set_guesses_for(m.fs.load_sx[4].aqueous_inlet, tear_guesses1)
# seq.set_guesses_for(m.fs.aq_inter_mixer[2].sx, tear_guesses2)
# seq.set_guesses_for(m.fs.aq_inter_mixer[3].sx, tear_guesses2)

sx_initializer = MixerSettlerExtractionInitializer()
load_sx_units = [m.fs.load_sx[s] for s in load_stage_list]
# scrub_sx_units = m.fs.scrub_sx
# strip_sx_units = [m.fs.strip_sx[s] for s in strip_stage_list]

mixer_initializer = MixerInitializer()
load_interstage_mixer = [m.fs.org_inter_mixer[s] for s in load_interstage_list]
# strip_interstage_mixer = [m.fs.aq_inter_mixer[s] for s in strip_stage_list]

tank_initializer = BlockTriangularizationInitializer()
neutral_tank = m.fs.aq_feed_neutral


def function(unit):
    if unit in load_sx_units:
        print(degrees_of_freedom(unit))
        print(f"Initializing {unit}")
        sx_initializer.initialize(unit)
    # elif unit in scrub_sx_units:
    #     print(f"Initializing {unit}")
    #     sx_initializer.initialize(unit)
    # elif unit in strip_sx_units:
    #     print(f"Initializing {unit}")
    #     sx_initializer.initialize(unit)
    elif unit in load_interstage_mixer:
        print(degrees_of_freedom(unit))
        print(f"Initializing {unit}")
        mixer_initializer.initialize(unit)
    # elif unit in strip_interstage_mixer:
    #     print(degrees_of_freedom(unit))
    #     print(f"Initializing {unit}")
    #     mixer_initializer.initialize(unit)
    else:
        print(f"Initializing {unit}")
        tank_initializer.initialize(unit)


seq.run(m, function)

# solver = get_solver("ipopt_v2")
# solver.options["max_iter"] = 2000
# # solver.options["halt_on_ampl_error"] = "yes"
# # solver.options["nlp_scaling_method"] = "user-scaling"
# solver.solve(m, tee=True)
