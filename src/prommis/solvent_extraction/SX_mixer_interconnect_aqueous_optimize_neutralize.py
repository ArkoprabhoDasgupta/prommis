from pyomo.environ import (
    ConcreteModel,
    units,
    RangeSet,
    TransformationFactory,
    Var,
    maximize,
    Param,
    value,
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

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution_new import REESolExOgParameters
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
    MixerSettlerExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package_new import (
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
    SingleControlVolumeUnitInitializer,
)

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()

dosage = 4

number_of_stages = 5
number_of_interstage_mixers = number_of_stages - 1

stage_list = RangeSet(1, number_of_stages)
interstage_list = RangeSet(1, number_of_interstage_mixers)

m.fs.mixer_settler_sx = MixerSettlerExtraction(
    stage_list,
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

m.fs.neutral = NeutralizationTank(interstage_list, property_package=m.fs.leach_soln)

# m.fs.interstage_mixer = Mixer(
#     interstage_list,
#     property_package=m.fs.leach_soln,
#     num_inlets=2,
#     inlet_list=["aqueous_out", "feed"],
#     material_balance_type=MaterialBalanceType.componentTotal,
#     energy_mixing_type=MixingType.none,
#     momentum_mixing_type=MomentumMixingType.none,
# )

pH = 1.6
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "H"].fix(
    10 ** (-pH) * 1e3
)  # decision variable
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "SO4"].fix(10 ** (-pH) * 96e3)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-5)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Al"].fix(422.375)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Ca"].fix(109.542)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Fe"].fix(688.266)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.032)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.124)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "La"].fix(0.986)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Ce"].fix(2.277)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.303)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Nd"].fix(0.946)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.097)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.2584)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.047)

m.fs.mixer_settler_sx[1].aqueous_inlet.flow_vol.fix(62.01)

m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Kerosene"].fix(
    820e3
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
    975e3 * 5 / 100
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Al_o"].fix(
    1.267e-5
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Ca_o"].fix(
    2.684e-5
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Fe_o"].fix(
    2.873e-6
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Sc_o"].fix(
    1.734
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Y_o"].fix(
    2.179e-5
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "La_o"].fix(
    0.000105
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Ce_o"].fix(
    0.00031
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Pr_o"].fix(
    3.711e-5
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Nd_o"].fix(
    0.000165
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Sm_o"].fix(
    1.701e-5
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Gd_o"].fix(
    3.357e-5
)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Dy_o"].fix(
    8.008e-6
)

m.fs.mixer_settler_sx[number_of_stages].organic_inlet.flow_vol.fix(62.01)

m.fs.mixer_settler_sx[:].mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.mixer_settler_sx[:].mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
    305.15 * units.K
)

m.fs.mixer_settler_sx[:].mixer[:].unit.mscontactor.volume[:].fix(1e-2 * units.m**3)

m.fs.mixer_settler_sx[:].organic_settler[:].unit.area.fix(1e-2)
m.fs.mixer_settler_sx[:].aqueous_settler[:].unit.area.fix(1e-2)
m.fs.mixer_settler_sx[:].aqueous_settler[:].unit.length.fix(1e-2)
m.fs.mixer_settler_sx[:].organic_settler[:].unit.length.fix(1e-2)


for i in stage_list:
    if i != 1:
        m.add_component(
            f"organic_sx_{i}_to_{i-1}",
            Arc(
                source=m.fs.mixer_settler_sx[i].organic_outlet,
                destination=m.fs.mixer_settler_sx[i - 1].organic_inlet,
            ),
        )

for i in interstage_list:
    m.add_component(
        f"neutral_tank_{i}_to_aqueous_sx_{i+1}",
        Arc(
            source=m.fs.neutral[i].outlet,
            destination=m.fs.mixer_settler_sx[i + 1].aqueous_inlet,
        ),
    )
    m.add_component(
        f"aqueous_sx_{i}_to_neutral_tank_{i}",
        Arc(
            source=m.fs.mixer_settler_sx[i].aqueous_outlet,
            destination=m.fs.neutral[i].inlet,
        ),
    )

TransformationFactory("network.expand_arcs").apply_to(m)


m.fs.neutral[:].base_flowrate[0].fix(1)
m.fs.neutral[:].base_concentration[0].fix(0.01)

m.fs.neutral[:].control_volume.properties_out[0.0].temperature.fix(305)
m.fs.neutral[:].control_volume.properties_out[0.0].pressure.fix(101325)


# # trial_element_list = ["Y", "Dy", "Gd", "Sm", "Nd", "Ce"]

# # m.percentage_recovery = Var(trial_element_list, initialize=1, bounds=(0, 100))

# m.tree_recovery = Var(initialize=1, bounds=(0, 100))

# # @m.Constraint(trial_element_list)
# # def recovery_constraint(m, e):
# #     return (
# #         m.percentage_recovery[e]
# #         == (
# #             (
# #                 m.fs.mixer_settler_sx[1].organic_outlet.conc_mass_comp[0, f"{e}_o"]
# #                 * m.fs.mixer_settler_sx[1].organic_outlet.flow_vol[0]
# #                 - m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[
# #                     0, f"{e}_o"
# #                 ]
# #                 * m.fs.mixer_settler_sx[number_of_stages].organic_inlet.flow_vol[0]
# #             )
# #             / (
# #                 m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, f"{e}"]
# #                 * m.fs.mixer_settler_sx[1].aqueous_inlet.flow_vol[0]
# #             )
# #         )
# #         * 100
# #     )


# @m.Constraint()
# def recovery_constraint(m):
#     return (
#         m.tree_recovery
#         == (
#             (
#                 sum(
#                     m.fs.mixer_settler_sx[1].organic_outlet.conc_mass_comp[0, f"{e}_o"]
#                     for e in m.fs.leach_soln.component_list
#                     if e
#                     not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
#                 )
#                 * m.fs.mixer_settler_sx[1].organic_outlet.flow_vol[0]
#                 - sum(
#                     m.fs.mixer_settler_sx[
#                         number_of_stages
#                     ].organic_inlet.conc_mass_comp[0, f"{e}_o"]
#                     for e in m.fs.leach_soln.component_list
#                     if e
#                     not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
#                 )
#                 * m.fs.mixer_settler_sx[number_of_stages].organic_inlet.flow_vol[0]
#             )
#             / (
#                 sum(
#                     m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, e]
#                     for e in m.fs.leach_soln.component_list
#                     if e
#                     not in ["H2O", "H", "HSO4", "SO4", "Cl", "Al", "Ca", "Fe", "Sc"]
#                 )
#                 * m.fs.mixer_settler_sx[1].aqueous_inlet.flow_vol[0]
#             )
#         )
#         * 100
#     )


# REE_set = ["Ce_o", "Dy_o", "Gd_o", "La_o", "Nd_o", "Pr_o", "Sm_o", "Y_o"]

# C_organic_max = {
#     "Ce_o": 0.205,
#     "Dy_o": 0.047,
#     "Gd_o": 0.125,
#     "La_o": 0.093,
#     "Nd_o": 0.068,
#     "Pr_o": 0.0307,
#     "Sm_o": 0.0196,
#     "Y_o": 0.124,
# }


# @m.Constraint()
# def sx_inlet_H_constraint(m):
#     return (
#         m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "H"]
#         >= 3 * units.mg / units.L
#     )


# @m.Constraint(interstage_list)
# def interstage_inlet_H_constraint(m, s):
#     return (
#         m.fs.interstage_mixer[s].feed.conc_mass_comp[0, "H"] >= 3 * units.mg / units.L
#     )


# @m.Constraint(REE_set)
# def REE_conc_constraint(m, e):
#     return (
#         m.fs.mixer_settler_sx[1].organic_outlet.conc_mass_comp[0, e] >= C_organic_max[e]
#     )


# # m.rho = Param(initialize=5)


# @m.Objective(sense=maximize)
# def objective_function(m):
#     return m.tree_recovery


print(degrees_of_freedom(m))

seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
# seq.options.iterLim = 3

# Using the SD tool
G = seq.create_graph(m)
heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
order = seq.calculation_order(G)

tear_guesses1 = {
    "flow_vol": {0: 62.01},
    "conc_mass_comp": {
        (0, "Al_o"): 0.048,
        (0, "Ca_o"): 1.98e-2,
        (0, "Ce_o"): 5.71e-3,
        (0, "Dy_o"): 1.077,
        (0, "Fe_o"): 1.954,
        (0, "Gd_o"): 0.14,
        (0, "La_o"): 4.03e-3,
        (0, "Nd_o"): 3.37e-3,
        (0, "Pr_o"): 1.04e-3,
        (0, "Sc_o"): 1.74,
        (0, "Sm_o"): 4.91e-3,
        (0, "Y_o"): 4.17,
        (0, "DEHPA"): 9.8e5 * 0.05,
        (0, "Kerosene"): 8.2e5,
    },
    "temperature": {0: 305},
    "pressure": {0: 101325},
}
tear_guesses2 = {
    "flow_vol": {0: 62.07},
    "conc_mass_comp": {
        (0, "Al"): 320.46,
        (0, "Ca"): 62.14,
        (0, "Ce"): 3.26,
        (0, "Cl"): 192.63,
        (0, "Dy"): 4.6e-2,
        (0, "Fe"): 452.28,
        (0, "Gd"): 0.40,
        (0, "H"): 2.92,
        (0, "H2O"): 1000000,
        (0, "HSO4"): 732.71,
        (0, "La"): 1.18,
        (0, "Nd"): 1.63,
        (0, "Pr"): 0.41,
        (0, "SO4"): 2543.95,
        (0, "Sc"): 2.25e-2,
        (0, "Sm"): 0.16,
        (0, "Y"): 0.11,
    },
    "temperature": {0: 305},
    "pressure": {0: 101325},
}

seq.set_guesses_for(m.fs.mixer_settler_sx[1].organic_inlet, tear_guesses1)
seq.set_guesses_for(m.fs.neutral[2].inlet, tear_guesses2)
seq.set_guesses_for(m.fs.neutral[3].inlet, tear_guesses2)
seq.set_guesses_for(m.fs.neutral[4].inlet, tear_guesses2)

mx_sx = [
    m.fs.mixer_settler_sx[1],
    m.fs.mixer_settler_sx[2],
    m.fs.mixer_settler_sx[3],
    m.fs.mixer_settler_sx[4],
    m.fs.mixer_settler_sx[5],
]

neutral_tank = [
    m.fs.neutral[1],
    m.fs.neutral[2],
    m.fs.neutral[3],
    m.fs.neutral[4],
]


def function(unit):
    if unit in mx_sx:
        MixerSettlerExtractionInitializer().initialize(unit)
        print(f"Initialized {unit}")
    elif unit in neutral_tank:
        BlockTriangularizationInitializer().initialize(unit)
        print(f"Initialized {unit}")


seq.run(m, function)

# solver = get_solver("ipopt_v2")
# solver.options["max_iter"] = 10000
# solver.solve(m, tee=True)
