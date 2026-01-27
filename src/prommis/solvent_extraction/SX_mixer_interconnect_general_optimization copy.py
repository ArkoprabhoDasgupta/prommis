from pyomo.environ import (
    ConcreteModel,
    units,
    RangeSet,
    TransformationFactory,
    Var,
    maximize,
    Param,
    value,
    minimize, Constraint, Suffix
)
from pyomo.network import Arc

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
)
from idaes.core.util import to_json
from idaes.core.util.model_diagnostics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import set_scaling_factor

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution_new import REESolExOgParameters
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
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

m.fs.aq_inter_mixer = Mixer(
    strip_stage_list,
    property_package=m.fs.leach_soln,
    num_inlets=2,
    inlet_list=["sx", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

# m.fs.aq_scrub_mixer = Mixer(
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
m.neutral_to_aq_sx = Arc(
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
    if i != loading_stages:
        m.add_component(
            f"load_interstage_organic_{i}_to_sx_organic_{i}",
            Arc(
                source=m.fs.org_inter_mixer[i].outlet,
                destination=m.fs.load_sx[i].organic_inlet,
            ),
        )

# stripping arcs

for i in strip_stage_list:
    # stripping aqueous mixer to aqueous sx
    m.add_component(
        f"strip_aqueous_mixer_{i}_to_sx_{i}",
        Arc(
            source=m.fs.aq_inter_mixer[i].outlet,
            destination=m.fs.strip_sx[i].aqueous_inlet,
        ),
    )
    if i != strip_stages:
        # stripping sx mixer to aqueous sx
        m.add_component(
            f"strip_aqueous_sx_{i}_to_mixer_{i+1}",
            Arc(
                source=m.fs.strip_sx[i].aqueous_outlet,
                destination=m.fs.aq_inter_mixer[i + 1].sx,
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


# connect scrub mixer and scrub sx

# m.aqueous_scrub_mixer_to_sx = Arc(
#     source=m.fs.aq_scrub_mixer.outlet,
#     destination=m.fs.scrub_sx.aqueous_inlet,
# )

# connect organic loading and scrubbing and stripping

m.organic_load_to_scrub = Arc(
    source=m.fs.load_sx[1].organic_outlet,
    destination=m.fs.scrub_sx.organic_inlet,
)

m.organic_scrub_to_strip = Arc(
    source=m.fs.scrub_sx.organic_outlet,
    destination=m.fs.strip_sx[strip_stages].organic_inlet,
)

TransformationFactory("network.expand_arcs").apply_to(m)

# set inlets

# neutralization tank input

pH_load = 1
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
# m.fs.aq_feed_neutral.base_flowrate[0].fix(1)
# m.fs.aq_feed_neutral.base_concentration[0].fix(0.01) #dv

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
    # dosage = 10  # in vol %
    # m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "DEHPA"].fix(975.8e3 * dosage / 100) # dv

    # m.fs.org_inter_mixer[i].feed.flow_vol.setub(5)
    # m.fs.org_inter_mixer[i].feed.flow_vol.setlb(1)

# stripping aqueous interstage addition

for i in strip_stage_list:
    for e in m.fs.leach_soln.component_list:
        if e not in ["H2O", "HSO4", "SO4", "H"]:
            m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-9)
    m.fs.aq_inter_mixer[i].feed.flow_vol.fix(3)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H2O"].fix(1e6)
    # m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].fix(
    #     10 ** (-pH_strip) * 2 * units.gram / units.L
    # ) #dv
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "SO4"].fix(1e-8)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "HSO4"].fix(1e-8)

# aqueous strip sx inlet

pH_strip = 0.4

for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "HSO4", "SO4", "H"]:
        m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, e].fix(1e-9)
m.fs.aq_inter_mixer[1].sx.flow_vol.fix(62.01)
m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "H"].fix(
    1 * units.gram / units.L
)  # maybe dv
m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "SO4"].fix(1e-8)
m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "HSO4"].fix(1e-8)

# # aqueous scrub sx inlet

# for e in m.fs.leach_soln.component_list:
#     if e not in ["H2O", "HSO4", "SO4", "H"]:
#         m.fs.aq_scrub_mixer.sx.conc_mass_comp[0, e].fix(1e-9)
# m.fs.aq_scrub_mixer.sx.flow_vol.fix(62.01)
# m.fs.aq_scrub_mixer.sx.conc_mass_comp[0, "H2O"].fix(1e6)
# m.fs.aq_scrub_mixer.sx.conc_mass_comp[0, "H"].fix(
#     0.1 * units.gram / units.L
# )  # maybe dv
# m.fs.aq_scrub_mixer.sx.conc_mass_comp[0, "SO4"].fix(1e-8)
# m.fs.aq_scrub_mixer.sx.conc_mass_comp[0, "HSO4"].fix(1e-8)

# aqueous scrub sx feed

for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "HSO4", "SO4", "H"]:
        m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, e].fix(1e-9)
m.fs.scrub_sx.aqueous_inlet.flow_vol.fix(1)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H"].fix(
    0.1 * units.gram / units.L
) #dv
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(1e-8)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-8)


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

m.fs.strip_sx[:].mixer[:].unit.mscontactor.volume[:].fix(1e-2 * units.m**3)
m.fs.strip_sx[:].organic_settler[:].unit.area.fix(1e-2)
m.fs.strip_sx[:].aqueous_settler[:].unit.area.fix(1e-2)
m.fs.strip_sx[:].aqueous_settler[:].unit.length.fix(1e-2)
m.fs.strip_sx[:].organic_settler[:].unit.length.fix(1e-2)
m.fs.strip_sx[:].mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.strip_sx[:].mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
    305.15 * units.K
)

m.fs.scrub_sx.mixer[:].unit.mscontactor.volume[:].fix(1e-2 * units.m**3)
m.fs.scrub_sx.organic_settler[:].unit.area.fix(1e-2)
m.fs.scrub_sx.aqueous_settler[:].unit.area.fix(1e-2)
m.fs.scrub_sx.aqueous_settler[:].unit.length.fix(1e-2)
m.fs.scrub_sx.organic_settler[:].unit.length.fix(1e-2)
m.fs.scrub_sx.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.scrub_sx.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
    305.15 * units.K
)

# fix neutral tank volumes and temperatures

m.fs.aq_feed_neutral.control_volume.properties_out[0.0].temperature.fix(305)
m.fs.aq_feed_neutral.control_volume.properties_out[0.0].pressure.fix(101325)

# fix interstage mixer

m.fs.org_inter_mixer[:].mixed_state[0.0].temperature.fix(305)
m.fs.aq_inter_mixer[:].mixed_state[0.0].temperature.fix(305)
m.fs.org_inter_mixer[:].mixed_state[0.0].pressure.fix(101325)
m.fs.aq_inter_mixer[:].mixed_state[0.0].pressure.fix(101325)
# m.fs.aq_scrub_mixer.mixed_state[0.0].pressure.fix(101325)
# m.fs.aq_scrub_mixer.mixed_state[0.0].temperature.fix(305)


m.fs.strip_sx[:].mixer[1].unit.mscontactor.heterogeneous_reactions[0.0,1].ascorbic_dosage.fix(1e-8)
m.fs.scrub_sx.mixer[1].unit.mscontactor.heterogeneous_reactions[0.0,1].ascorbic_dosage.fix(1e-8)

@m.Constraint(load_stage_list)
def ascorbic_acid_constraint(m, s):
    if s==1:
        return Constraint.Skip
    else:
        return m.fs.load_sx[s].mixer[1].unit.mscontactor.heterogeneous_reactions[0.0,1].ascorbic_dosage == m.fs.load_sx[s-1].mixer[1].unit.mscontactor.heterogeneous_reactions[0.0,1].ascorbic_dosage

@m.Constraint(load_stage_list)
def organic_dosage_constraint(m, s):
    if s==loading_stages:
        return Constraint.Skip
    else:
        return m.fs.load_sx[s].organic_inlet.conc_mass_comp[0, "DEHPA"]/9758 >= m.fs.load_sx[s+1].organic_inlet.conc_mass_comp[0, "DEHPA"]/9758


print(degrees_of_freedom(m))
# assert 1==2

REE_list = [
    e
    for e in m.fs.leach_soln.component_list
    if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Al", "Fe", "Ca"]
]

m.percentage_recovery = Var(REE_list, initialize=0.9, bounds=(0, 1))


@m.Constraint(REE_list)
def recovery_constraint(m, e):
    scrubbing_inlet = (
        m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, e] * m.fs.scrub_sx.aqueous_inlet.flow_vol[0]
    )
    scrubbing_outlet = (
        m.fs.scrub_sx.aqueous_outlet.conc_mass_comp[0, e] * m.fs.scrub_sx.aqueous_outlet.flow_vol[0]
    )
    stripping_inlet = (
        m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, e] * m.fs.aq_inter_mixer[1].sx.flow_vol[0]
    )   
    stripping_outlet = (
        m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, e] * m.fs.strip_sx[strip_stages].aqueous_outlet.flow_vol[0]
    )
    load_inlet = (
        m.fs.aq_feed_neutral.inlet.flow_vol[0] * m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, e]
    )   
    return m.percentage_recovery[e] == (scrubbing_outlet - scrubbing_inlet + stripping_outlet - stripping_inlet) / load_inlet


m.tree_recovery = Var(initialize=0.3, bounds=(0, 1))

@m.Constraint()
def tree_recovery_constraint(m):

    scrubbing_inlet = sum(m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, e] for e in REE_list) * m.fs.scrub_sx.aqueous_inlet.flow_vol[0]
    scrubbing_outlet = sum(m.fs.scrub_sx.aqueous_outlet.conc_mass_comp[0, e] for e in REE_list) * m.fs.scrub_sx.aqueous_outlet.flow_vol[0]

    stripping_inlet = sum(m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, e] for e in REE_list) * m.fs.aq_inter_mixer[1].sx.flow_vol[0]
    stripping_outlet = sum(m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, e] for e in REE_list) * m.fs.strip_sx[strip_stages].aqueous_outlet.flow_vol[0]

    load_inlet = m.fs.aq_feed_neutral.inlet.flow_vol[0] * sum(m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, e] for e in REE_list)
    
    
    return m.tree_recovery == (scrubbing_outlet - scrubbing_inlet + stripping_outlet - stripping_inlet)/ load_inlet

# @m.Constraint()
# def tree_recovery_constraint(m):
#     load_inlet = m.fs.aq_feed_neutral.inlet.flow_vol[0] * sum(m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, e] for e in REE_list)
#     load_outlet = m.fs.load_sx[loading_stages].aqueous_outlet.flow_vol[0] * sum(m.fs.load_sx[loading_stages].aqueous_outlet.conc_mass_comp[0, e] for e in REE_list)
#     return m.tree_recovery == (1-load_outlet/load_inlet)*100

# @m.Objective(sense=maximize)
# def objective_function(m):
#     return m.tree_recovery


# set upper bounds to decision variables
m.fs.aq_feed_neutral.base_concentration.setlb(0.01)
# m.fs.aq_feed_neutral.base_concentration.fix(3)
m.fs.aq_feed_neutral.base_flowrate.setub(5)
m.fs.aq_feed_neutral.base_flowrate.setlb(1)

for i in strip_stage_list:
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].setub(2 * units.gram / units.L)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].setlb(
        1e-1 * units.gram / units.L
    )
# m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "H"].setub(2 * units.gram / units.L)

m.product_distribution = Var(REE_list, initialize=0.12, bounds=(0, 1))

# @m.Constraint(m.fs.time)
# def neutral_base_restriction(m,t):
#     return


@m.Constraint(REE_list)
def product_distribution_constraint(m, e):
    return m.product_distribution[e] == (
        m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, e]
        / sum(
            m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, s]
            for s in REE_list
        )
    )


product_distribution = {
    "Y": 23,
    "La": 9,
    "Ce": 26,
    "Pr": 3,
    "Nd": 16,
    "Sm": 5,
    "Gd": 6,
    "Dy": 5,
}

# market_demand = {
#     "Y": 19,
#     "La": 8,
#     "Ce": 26,
#     "Pr": 4,
#     "Nd": 20,
#     "Sm": 6,
#     "Gd": 6,
#     "Dy": 4,
# }


m.R_I_dist = Var(["ree", "Al", "Ca", "Fe", "Sc"], bounds=(0, 1), initialize=0.5)

@m.Constraint()
def ree_composition(m):
    return m.R_I_dist["ree"] == sum(
        m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, e]
        for e in REE_list
    ) / sum(
        m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, e]
        for e in m.fs.leach_soln.component_list
        if e not in ["H2O", "H", "SO4", "HSO4", "Cl"]
    )


@m.Constraint(["Al", "Ca", "Fe", "Sc"])
def impurity_composition(m, e):
    return m.R_I_dist[e] == m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[
        0, e
    ] / sum(
        m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, e]
        for e in m.fs.leach_soln.component_list
        if e not in ["H2O", "H", "SO4", "HSO4", "Cl"]
    )


# @m.Constraint()
# def production_constraint(m):
#     return m.R_I_dist["Al"] + m.R_I_dist["Fe"] <= 0.6


# m.Y_Ce_combo = Var(initialize=0.5, bounds=(0, 1))


# @m.Constraint()
# def Y_Ce_constraint(m):
#     return m.Y_Ce_combo == (
#         sum(
#             m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, e]
#             for e in ["Y", "Ce"]
#         )
#         * m.fs.strip_sx[strip_stages].aqueous_outlet.flow_vol[0]
#         - sum(m.fs.aq_inter_mixer[2].sx.conc_mass_comp[0, e] for e in ["Y", "Ce"])
#         * m.fs.aq_inter_mixer[2].sx.flow_vol[0]
#     ) / (
#         m.fs.aq_feed_neutral.inlet.flow_vol[0]
#         * sum(m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, e] for e in ["Y", "Ce"])
#     )


# m.percentage_recovery["Gd"].setub(30)
# m.R_I_dist["Al"].setub(0.6)



@m.Objective(sense=maximize)
def objective_function(m):
    return m.tree_recovery
    # return m.tree_recovery + 1e-3 * (0.6 - m.R_I_dist["Al"] - m.R_I_dist["Fe"])
    # return m.percentage_recovery["Gd"] + 0.03 * (
    #     10 - m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, "Y"]
    # )
    # return m.percentage_recovery["Gd"] + 40 * (0.12 - m.product_distribution["Ce"])
    # return m.percentage_recovery["Y"] + 1.2e-2 * (0.13 - m.product_distribution["Ce"])
    # return m.Y_Ce_combo
    # return m.R_I_dist["Al"] + m.R_I_dist["Fe"]

m.scaling_factor = Suffix(direction=Suffix.EXPORT)

# set_scaling_factor(m.fs.aq_inter_mixer[1].feed_state[0.0].pH_constraint['liquid'], 1e3)
for s in strip_stage_list:
    if s!=2:
        set_scaling_factor(m.fs.aq_inter_mixer[s].feed_state[0.0].pH_constraint['liquid'], 1e-2)
    else:
        set_scaling_factor(m.fs.aq_inter_mixer[s].sx_state[0.0].pH_constraint['liquid'], 1e3)
    for e in REE_list:
        set_scaling_factor(m.fs.strip_sx[s].mixer[1].unit.distribution_extent_constraint[0,1,e], 1)
    set_scaling_factor(m.fs.strip_sx[s].mixer[1].unit.distribution_extent_constraint[0,1,'Fe'], 1)
    set_scaling_factor(m.fs.strip_sx[s].mixer[1].unit.distribution_extent_constraint[0,1,'Al'], 1)
for s in load_stage_list:
    for e in REE_list:
        set_scaling_factor(m.fs.load_sx[s].mixer[1].unit.distribution_extent_constraint[0,1,e], 1)


scaling = TransformationFactory("core.scale_model")
scaled_model = scaling.create_using(m, rename=False)

print(degrees_of_freedom(scaled_model))

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 20000
# solver.options["halt_on_ampl_error"] = "yes"
# solver.options["nlp_scaling_method"] = "user-scaling"
solver.solve(scaled_model, tee=True)

scaling.propagate_solution(scaled_model, m)

decision_vars = {
    "Neutralization tank base flowrate": value(m.fs.aq_feed_neutral.base_flowrate[0]),
    "Neutralization tank base concentration": value(m.fs.aq_feed_neutral.base_concentration[0]),
    'Organic interstage mixer dosages': [value(m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "DEHPA"])/9758 for i in range(1, loading_stages)],
    "Aqueous interstage mixer acidities": [value(m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"]) for i in range(1, strip_stages + 1)],
}