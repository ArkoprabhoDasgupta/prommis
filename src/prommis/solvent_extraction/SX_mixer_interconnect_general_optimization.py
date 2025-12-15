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

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()

loading_stages = 4

load_stage_list = RangeSet(1, loading_stages)
load_interstage_list = RangeSet(1, loading_stages - 1)

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

m.fs.org_inter_mixer = Mixer(
    load_interstage_list,
    property_package=m.fs.prop_o,
    num_inlets=2,
    inlet_list=["sx", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

m.fs.aq_feed_neutral = NeutralizationTank(property_package=m.fs.leach_soln)

strip_stages = 4
strip_stage_list = RangeSet(1, strip_stages)

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

m.fs.aq_inter_mixer = Mixer(
    strip_stage_list,
    property_package=m.fs.leach_soln,
    num_inlets=2,
    inlet_list=["sx", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

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

# connect organic loading and stripping
m.organic_load_to_strip = Arc(
    source=m.fs.load_sx[1].organic_outlet,
    destination=m.fs.strip_sx[strip_stages].organic_inlet,
)

TransformationFactory("network.expand_arcs").apply_to(m)

pH_load = 1
m.fs.aq_feed_neutral.inlet.flow_vol.fix(62.01)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Al"].fix(400)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ca"].fix(100)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Fe"].fix(600)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Sc"].fix(16.25)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "La"].fix(62.37)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ce"].fix(132.5)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Pr"].fix(15.36)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Nd"].fix(59.19)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Sm"].fix(10.29)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Gd"].fix(8.44)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Dy"].fix(6.84)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Y"].fix(33.45)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "H"].fix(
    10 ** (-pH_load) * units.gram / units.L
)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "SO4"].fix(
    10 ** (-pH_load) * 48 * units.gram / units.L
)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "HSO4"].fix(1e-4)

for e in m.fs.prop_o.component_list:
    if e not in ["Kerosene", "DEHPA"]:
        m.fs.load_sx[loading_stages].organic_inlet.conc_mass_comp[0, e].fix(1e-9)
m.fs.load_sx[loading_stages].organic_inlet.flow_vol.fix(62.01)
m.fs.load_sx[loading_stages].organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)
dosage = 8  # in vol %
m.fs.load_sx[loading_stages].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
    975.8e3 * dosage / 100
)

for i in load_interstage_list:
    for e in m.fs.prop_o.component_list:
        if e not in ["Kerosene", "DEHPA"]:
            m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-9)
    m.fs.org_inter_mixer[i].feed.flow_vol.fix(1)
    m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "Kerosene"].fix(820e3)
    # dosage = 10  # in vol %
    # m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "DEHPA"].fix(975.8e3 * dosage / 100) # dv

pH_strip = 0.5

for i in strip_stage_list:
    for e in m.fs.leach_soln.component_list:
        if e not in ["H2O", "HSO4", "SO4", "H"]:
            m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-9)
    m.fs.aq_inter_mixer[i].feed.flow_vol.fix(1)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H2O"].fix(1e6)
    # m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].fix(
    #     10 ** (-pH_strip) * 2 * units.gram / units.L
    # ) #dv
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "SO4"].fix(1e-8)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "HSO4"].fix(1e-8)

# aqueous strip inlet
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "HSO4", "SO4", "H"]:
        m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, e].fix(1e-9)
m.fs.aq_inter_mixer[1].sx.flow_vol.fix(62.01)
m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "H"].fix(
    10 ** (-pH_strip) * 1 * units.gram / units.L
)  # maybe dv
m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "SO4"].fix(1e-8)
m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "HSO4"].fix(1e-8)


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

# fix neutral tank volumes and temperatures
m.fs.aq_feed_neutral.base_flowrate[0].fix(2)
# m.fs.aq_feed_neutral.base_concentration[0].fix(0.01) #dv
m.fs.aq_feed_neutral.control_volume.properties_out[0.0].temperature.fix(305)
m.fs.aq_feed_neutral.control_volume.properties_out[0.0].pressure.fix(101325)

# fix interstage mixer
m.fs.org_inter_mixer[:].mixed_state[0.0].temperature.fix(305)
m.fs.aq_inter_mixer[:].mixed_state[0.0].temperature.fix(305)
m.fs.org_inter_mixer[:].mixed_state[0.0].pressure.fix(101325)
m.fs.aq_inter_mixer[:].mixed_state[0.0].pressure.fix(101325)

REE_list = [
    e
    for e in m.fs.leach_soln.component_list
    if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Al", "Fe", "Ca", "Sc"]
]

m.percentage_recovery = Var(REE_list, initialize=0.9, bounds=(0, 1))


@m.Constraint(REE_list)
def recovery_constraint(m, e):
    return m.percentage_recovery[e] == (
        (
            m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, e]
            * m.fs.strip_sx[strip_stages].aqueous_outlet.flow_vol[0]
            - m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, e]
            * m.fs.aq_inter_mixer[1].sx.flow_vol[0]
        )
        / (
            m.fs.aq_feed_neutral.inlet.flow_vol[0]
            * m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, e]
        )
    )


m.tree_recovery = Var(initialize=1, bounds=(0, 100))


@m.Constraint()
def tree_recovery_constraint(m):
    return m.tree_recovery == (
        (
            sum(
                m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, e]
                for e in REE_list
            )
            * m.fs.strip_sx[strip_stages].aqueous_outlet.flow_vol[0]
            - sum(m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, e] for e in REE_list)
            * m.fs.aq_inter_mixer[1].sx.flow_vol[0]
        )
        / (
            m.fs.aq_feed_neutral.inlet.flow_vol[0]
            * sum(m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, e] for e in REE_list)
        )
    )


# @m.Objective(sense=maximize)
# def objective_function(m):
#     return m.tree_recovery


# set upper bounds to decision variables
m.fs.aq_feed_neutral.base_concentration.setub(2)

for i in strip_stage_list:
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].setub(2 * units.gram / units.L)
    # m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].setlb(
    #     1e-4 * units.gram / units.L
    # )
m.fs.aq_inter_mixer[1].sx.conc_mass_comp[0, "H"].setub(2 * units.gram / units.L)

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

# product_constraint_element = ["Ce"]


# @m.Constraint()
# def production_constraint(m):
#     return m.product_distribution["Y"] <= 0.4

# m.product_distribution["Ce"].setub(0.13)

# @m.Constraint()
# def production_constraint(m):
#     return m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, "Ce"] <= 6


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


# m.percentage_recovery["Gd"].setub(30)
m.R_I_dist["Fe"].setub(0.1)


@m.Objective(sense=maximize)
def objective_function(m):
    # return m.tree_recovery
    return m.tree_recovery + 50 * (0.1 - m.R_I_dist["Fe"])
    # return m.percentage_recovery["Gd"] + 0.03 * (
    #     10 - m.fs.strip_sx[strip_stages].aqueous_outlet.conc_mass_comp[0, "Y"]
    # )
    # return m.percentage_recovery["Gd"] + 40 * (0.12 - m.product_distribution["Ce"])
    # return m.percentage_recovery["Y"] + 1.2e-2 * (0.13 - m.product_distribution["Ce"])
    # return m.percentage_recovery["Y"]


print(degrees_of_freedom(m))

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 30000
# solver.options["halt_on_ampl_error"] = "yes"
solver.solve(m, tee=True)
