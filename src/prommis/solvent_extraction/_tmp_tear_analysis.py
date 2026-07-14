from pyomo.environ import ConcreteModel, RangeSet, TransformationFactory
from idaes.core import FlowDirection, FlowsheetBlock, MaterialBalanceType
from pyomo.network import Arc, SequentialDecomposition
from prommis.solvent_extraction.leach_solution_properties_optimization import (
    LeachSolutionParameters,
)
from prommis.solvent_extraction.ree_og_distribution_new_optimization import (
    REESolExOgParameters,
)
from prommis.solvent_extraction.mixer_settler_extraction import MixerSettlerExtraction
from prommis.solvent_extraction.solvent_extraction_reaction_package_new_modified_for_flowsheet_optimization import (
    SolventExtractionReactions,
)
from prommis.solvent_extraction.neutralization_tank import NeutralizationTank
from idaes.models.unit_models.mixer import Mixer, MixingType, MomentumMixingType

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()

dosage = 8
load_number_of_stages = 4
load_stage_list = RangeSet(1, load_number_of_stages)
load_interstage_list = RangeSet(1, load_number_of_stages - 1)
strip_number_of_stages = 3
strip_stage_list = RangeSet(1, strip_number_of_stages)
strip_interstage_list = RangeSet(1, strip_number_of_stages - 1)
m.fs.reaxn.extractant_dosage = dosage

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
    settler_finite_elements=4,
)
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
    settler_finite_elements=4,
)
m.fs.aq_feed_neutral = NeutralizationTank(property_package=m.fs.leach_soln)
m.fs.org_inter_mixer = Mixer(
    load_interstage_list,
    property_package=m.fs.prop_o,
    num_inlets=2,
    inlet_list=["sx", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)
m.fs.aq_inter_mixer = Mixer(
    strip_interstage_list,
    property_package=m.fs.leach_soln,
    num_inlets=2,
    inlet_list=["sx", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

m.fs.neutral_to_sx = Arc(
    source=m.fs.aq_feed_neutral.outlet, destination=m.fs.load_sx[1].aqueous_inlet
)
for i in load_stage_list:
    if i != 1:
        m.add_component(
            f"load_aqueous_sx_{i-1}_to_{i}",
            Arc(
                source=m.fs.load_sx[i - 1].aqueous_outlet,
                destination=m.fs.load_sx[i].aqueous_inlet,
            ),
        )
        m.add_component(
            f"load_organic_sx_{i}_to_interstage_{i-1}",
            Arc(
                source=m.fs.load_sx[i].organic_outlet,
                destination=m.fs.org_inter_mixer[i - 1].sx,
            ),
        )
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
        m.add_component(
            f"strip_aqueous_sx_{i}_to_mixer_{i}",
            Arc(
                source=m.fs.strip_sx[i].aqueous_outlet,
                destination=m.fs.aq_inter_mixer[i].sx,
            ),
        )
        m.add_component(
            f"strip_aqueous_mixer_{i}_to_sx_{i+1}",
            Arc(
                source=m.fs.aq_inter_mixer[i].outlet,
                destination=m.fs.strip_sx[i + 1].aqueous_inlet,
            ),
        )
        m.add_component(
            f"strip_organic_sx_{i+1}_to_{i}",
            Arc(
                source=m.fs.strip_sx[i + 1].organic_outlet,
                destination=m.fs.strip_sx[i].organic_inlet,
            ),
        )
m.organic_load_to_scrub = Arc(
    source=m.fs.load_sx[1].organic_outlet, destination=m.fs.scrub_sx.organic_inlet
)
m.organic_scrub_to_strip = Arc(
    source=m.fs.scrub_sx.organic_outlet,
    destination=m.fs.strip_sx[strip_number_of_stages].organic_inlet,
)
TransformationFactory("network.expand_arcs").apply_to(m)

seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 3
G = seq.create_graph(m)
heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
print("tear_arcs:")
for a in heuristic_tear_set:
    print(a.name)
print("calculation_order:")
for n in seq.calculation_order(G):
    print(n)
