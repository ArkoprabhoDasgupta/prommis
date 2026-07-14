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
from pyomo.network import Arc
from prommis.solvent_extraction.leach_solution_properties_optimization import (
    LeachSolutionParameters,
)
from prommis.solvent_extraction.ree_og_distribution_new_optimization import (
    REESolExOgParameters,
)
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
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


m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()


# define sx

dosage = 8
number_of_stages = 4
stage_list = RangeSet(1, number_of_stages)
interstage_list = RangeSet(1, number_of_stages - 1)

m.fs.reaxn.extractant_dosage = dosage

m.fs.mixer_settler_ex = MixerSettlerExtraction(
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
    settler_finite_elements=2,
)

# define neutralization tank

m.fs.aq_feed_neutral = NeutralizationTank(property_package=m.fs.leach_soln)

# define ascorbic acid

m.fs.ascorbic_mixer = Mixer(
    property_package=m.fs.leach_soln,
    num_inlets=2,
    inlet_list=["ascorbic", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

m.fs.org_inter_mixer = Mixer(
    interstage_list,
    property_package=m.fs.prop_o,
    num_inlets=2,
    inlet_list=["sx", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

# define arcs

m.fs.neutral_to_asc_mixer = Arc(
    source=m.fs.aq_feed_neutral.outlet, destination=m.fs.ascorbic_mixer.feed
)

m.fs.asc_mixer_to_sx = Arc(
    source=m.fs.ascorbic_mixer.outlet,
    destination=m.fs.mixer_settler_ex[1].aqueous_inlet,
)

for i in stage_list:
    if i != 1:
        # loading sx aqueous phases
        m.add_component(
            f"load_aqueous_sx_{i-1}_to_{i}",
            Arc(
                source=m.fs.mixer_settler_ex[i - 1].aqueous_outlet,
                destination=m.fs.mixer_settler_ex[i].aqueous_inlet,
            ),
        )
        # loading sx organic to interstage
        m.add_component(
            f"load_organic_sx_{i}_to_interstage_{i-1}",
            Arc(
                source=m.fs.mixer_settler_ex[i].organic_outlet,
                destination=m.fs.org_inter_mixer[i - 1].sx,
            ),
        )
    # interstage to loading organic
    if i != number_of_stages:
        m.add_component(
            f"load_interstage_organic_{i}_to_sx_organic_{i}",
            Arc(
                source=m.fs.org_inter_mixer[i].outlet,
                destination=m.fs.mixer_settler_ex[i].organic_inlet,
            ),
        )

TransformationFactory("network.expand_arcs").apply_to(m)

# define inlet conditions

pH = 1.3
m.fs.aq_feed_neutral.inlet.flow_vol.fix(60.01)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ascorbic"].fix(1e-7)
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

m.fs.mixer_settler_ex[number_of_stages].organic_inlet.conc_mass_comp[0, "Kerosene"].fix(
    820e3
)
m.fs.mixer_settler_ex[number_of_stages].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
    975.8e3 * dosage / 100
)
for e in m.fs.prop_o.component_list:
    if e not in ["Kerosene", "DEHPA"]:
        m.fs.mixer_settler_ex[number_of_stages].organic_inlet.conc_mass_comp[0, e].fix(
            1e-7
        )
m.fs.mixer_settler_ex[number_of_stages].organic_inlet.flow_vol.fix(62.01)

m.fs.mixer_settler_ex[:].mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.mixer_settler_ex[:].mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
    305.15 * units.K
)
m.fs.mixer_settler_ex[:].mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.mixer_settler_ex[:].organic_settler[:].unit.area.fix(1)
m.fs.mixer_settler_ex[:].aqueous_settler[:].unit.area.fix(1)
m.fs.mixer_settler_ex[:].aqueous_settler[:].unit.length.fix(1)
m.fs.mixer_settler_ex[:].organic_settler[:].unit.length.fix(1)


m.fs.ascorbic_mixer.ascorbic.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.ascorbic_mixer.ascorbic.conc_mass_comp[0, "Ascorbic"].fix(40 * units.g / units.L)
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "Ascorbic"]:
        m.fs.ascorbic_mixer.ascorbic.conc_mass_comp[0, e].fix(1e-7)
m.fs.ascorbic_mixer.ascorbic.flow_vol.fix(4)

m.fs.ascorbic_mixer.mixed_state[0.0].pressure.fix(101235)
m.fs.ascorbic_mixer.mixed_state[0.0].temperature.fix(303.5)
m.fs.aq_feed_neutral.control_volume.properties_out[0.0].pressure.fix(101235)
m.fs.aq_feed_neutral.control_volume.properties_out[0.0].temperature.fix(303.5)

for i in interstage_list:
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


print(dof(m))

# # initializer = m.fs.mixer_settler_ex.default_initializer()
# # initializer.initialize(m.fs.mixer_settler_ex)

solver = get_solver("ipopt_v2")
solver.options["halt_on_ampl_error"] = "yes"
results = solver.solve(m, tee=True)

# percentage_recovery = {}
# for e in m.fs.leach_soln.component_list:
#     if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic"]:
#         percentage_recovery[e] = [
#             (
#                 (
#                     m.fs.mixer_settler_ex[1].organic_outlet.conc_mass_comp[
#                         t, f"{e}_o"
#                     ]()
#                     * m.fs.mixer_settler_ex[1].organic_outlet.flow_vol[t]()
#                     - m.fs.mixer_settler_ex[
#                         number_of_stages
#                     ].organic_inlet.conc_mass_comp[t, f"{e}_o"]()
#                     * m.fs.mixer_settler_ex[number_of_stages].organic_inlet.flow_vol[
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
#                 m.fs.mixer_settler_ex[1].organic_outlet.conc_mass_comp[t, f"{e}_o"]()
#                 * m.fs.mixer_settler_ex[1].organic_outlet.flow_vol[t]()
#                 for e in m.fs.leach_soln.component_list
#                 if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
#             )
#             - sum(
#                 m.fs.mixer_settler_ex[number_of_stages].organic_inlet.conc_mass_comp[
#                     t, f"{e}_o"
#                 ]()
#                 * m.fs.mixer_settler_ex[number_of_stages].organic_inlet.flow_vol[t]()
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
