#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import (
    ConcreteModel,
    maximize,
    units,
    TransformationFactory,
    RangeSet,
    value,
    Var,
    Constraint,
)

import os

import pandas as pd

from idaes.core import FlowDirection, FlowsheetBlock, MaterialBalanceType
from idaes.core.util import StoreSpec, to_json, from_json

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
from prommis.solvent_extraction.solvent_extraction import (
    SolventExtraction,
    SolventExtractionInitializer,
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
    SingleControlVolumeUnitInitializer,
)


m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()


# define stages
dosage = 16
load_number_of_stages = 3
load_stage_list = RangeSet(1, load_number_of_stages)
load_interstage_list = RangeSet(1, load_number_of_stages - 1)

strip_number_of_stages = 2
strip_stage_list = RangeSet(1, strip_number_of_stages)
strip_interstage_list = RangeSet(1, strip_number_of_stages - 1)

m.fs.reaxn.extractant_dosage = dosage

# define load sx
m.fs.load_sx = SolventExtraction(
    load_stage_list,
    number_of_finite_elements=1,
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
)

# define scrub sx
m.fs.scrub_sx = SolventExtraction(
    number_of_finite_elements=1,
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
)

# define strip sx
m.fs.strip_sx = SolventExtraction(
    strip_stage_list,
    number_of_finite_elements=1,
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
# m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ascorbic"].fix(2.44 * units.g / units.L)   df
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
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Cl"].fix(1e-6)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "H"].fix(10 ** (-pH) * units.g / units.L)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "SO4"].fix(
    10 ** (-pH) * 48 * units.g / units.L
)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "HSO4"].fix(1e-4)
# m.fs.aq_feed_neutral.base_flowrate[0].fix(0.5)     df
m.fs.aq_feed_neutral.base_concentration[0].fix(5)

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

m.fs.load_sx[:].mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
m.fs.load_sx[:].mscontactor.organic[:, :].temperature.fix(305.15 * units.K)
m.fs.load_sx[:].mscontactor.volume[:].fix(0.4 * units.m**3)


# define organic interstage inlet
for i in load_interstage_list:
    for e in m.fs.prop_o.component_list:
        if e not in ["Kerosene", "DEHPA"]:
            m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-7)
    # m.fs.org_inter_mixer[i].feed.flow_vol.fix(1)      df
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
            m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-7)
    # m.fs.aq_inter_mixer[i].feed.flow_vol.fix(0.5)             df
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H2O"].fix(1e6)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].fix(6 * units.g / units.L)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "Cl"].fix(
        6 * 35.5 * units.g / units.L
    )

m.fs.aq_inter_mixer[:].mixed_state[0.0].pressure.fix(101235)
m.fs.aq_inter_mixer[:].mixed_state[0.0].temperature.fix(303.5)

# define scrub sx aqueous inlet
# m.fs.scrub_sx.aqueous_inlet.flow_vol.fix(40.01)         df
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "Cl"]:
        m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, e].fix(1e-7)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H"].fix(0.1 * units.g / units.L)  # df
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(0.1 * 35.5 * units.g / units.L)

m.fs.scrub_sx.mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.scrub_sx.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
m.fs.scrub_sx.mscontactor.organic[:, :].temperature.fix(305.15 * units.K)


# define strip sx aqueous inlet
# m.fs.strip_sx[1].aqueous_inlet.flow_vol.fix(40.01)  # df
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "Cl"]:
        m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, e].fix(1e-7)
m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "H"].fix(6 * units.g / units.L)  # df
m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "Cl"].fix(6 * 35.5 * units.g / units.L)
m.fs.strip_sx[:].mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.strip_sx[:].mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
m.fs.strip_sx[:].mscontactor.organic[:, :].temperature.fix(305.15 * units.K)


@m.Constraint(load_stage_list)
def organic_dosage_constraint(m, s):
    if s == load_number_of_stages:
        return Constraint.Skip
    else:
        return (
            m.fs.load_sx[s].mscontactor.organic[0, 1].extractant_dosage
            <= m.fs.load_sx[s + 1].mscontactor.organic[0, 1].extractant_dosage
        )


@m.Constraint()
def sx_feed_pH_constraint(m):
    return (
        m.fs.aq_feed_neutral.control_volume.properties_out[0.0].pH_phase["liquid"] <= 3
    )


solver = get_solver("ipopt_v2")
solver.options["halt_on_ampl_error"] = "yes"
solver.options["max_iter"] = 4000

# m.Nd_Ce_recovery = Var(initialize=20, bounds=(0, 100))


# @m.Constraint()
# def Nd_Ce_recovery_constraint(m):

#     strip_inlet = sum(
#         m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, e]
#         * m.fs.strip_sx[1].aqueous_inlet.flow_vol[0]
#         for e in ["Ce", "Nd"]
#     )

#     scrub_inlet = sum(
#         m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, e]
#         * m.fs.scrub_sx.aqueous_inlet.flow_vol[0]
#         for e in ["Ce", "Nd"]
#     )

#     strip_outlet = sum(
#         m.fs.strip_sx[strip_number_of_stages].aqueous_outlet.conc_mass_comp[0, e]
#         * m.fs.strip_sx[strip_number_of_stages].aqueous_outlet.flow_vol[0]
#         for e in ["Ce", "Nd"]
#     )

#     scrub_outlet = sum(
#         m.fs.scrub_sx.aqueous_outlet.conc_mass_comp[0, e]
#         * m.fs.scrub_sx.aqueous_outlet.flow_vol[0]
#         for e in ["Ce", "Nd"]
#     )

#     feed_inlet = sum(
#         m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, e]
#         * m.fs.aq_feed_neutral.inlet.flow_vol[0]
#         for e in ["Ce", "Nd"]
#     )

#     # return (
#     #     100 * (strip_outlet + scrub_outlet - strip_inlet - scrub_inlet)
#     #     == m.Nd_Ce_recovery * feed_inlet
#     # )

#     return 100 * (strip_outlet - strip_inlet) == m.Nd_Ce_recovery * feed_inlet


m.tree_recovery = Var(initialize=0.3, bounds=(0, 1))


@m.Constraint()
def tree_recovery_constraint(m):

    strip_inlet = sum(
        m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, e]
        * m.fs.strip_sx[1].aqueous_inlet.flow_vol[0]
        for e in m.fs.leach_soln.component_list
        if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
    )

    strip_outlet = sum(
        m.fs.strip_sx[strip_number_of_stages].aqueous_outlet.conc_mass_comp[0, e]
        * m.fs.strip_sx[strip_number_of_stages].aqueous_outlet.flow_vol[0]
        for e in m.fs.leach_soln.component_list
        if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
    )

    feed_inlet = sum(
        m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, e]
        * m.fs.aq_feed_neutral.inlet.flow_vol[0]
        for e in m.fs.leach_soln.component_list
        if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
    )

    # return (
    #     100 * (strip_outlet + scrub_outlet - strip_inlet - scrub_inlet)
    #     == m.Nd_Ce_recovery * feed_inlet
    # )

    return (strip_outlet - strip_inlet) == m.tree_recovery * feed_inlet


# @m.Constraint()
# def Nd_Ce_composition_demand(m):

#     return m.Nd_Ce_composition >= 0.5


@m.Constraint()
def Fe_concentration_constraint(m):

    return (
        m.fs.strip_sx[strip_number_of_stages].aqueous_outlet.conc_mass_comp[0, "Fe"]
        <= 2.5e-4
    )


print(dof(m))

from_json(
    m,
    fname="sx_ascorbic_extraction_circuit_2.json",
    wts=StoreSpec.value(only_not_fixed=True),
)

# set bounds to the dfs
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ascorbic"].setlb(1.2 * units.g / units.L)
m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ascorbic"].setub(5 * units.g / units.L)

m.fs.aq_feed_neutral.base_flowrate[0].setlb(0.1)
m.fs.aq_feed_neutral.base_flowrate[0].setub(1)

for i in load_interstage_list:
    m.fs.org_inter_mixer[i].feed.flow_vol[0].setlb(0.1)
    m.fs.org_inter_mixer[i].feed.flow_vol[0].setub(3)

    # m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "DEHPA"].setlb(975.8e3 * 6 / 100)
    # m.fs.org_inter_mixer[i].feed.conc_mass_comp[0, "DEHPA"].setlb(975.8e3 * 20 / 100)


m.fs.scrub_sx.aqueous_inlet.flow_vol[0].setlb(20.01)
m.fs.scrub_sx.aqueous_inlet.flow_vol[0].setub(80.01)

m.fs.strip_sx[1].aqueous_inlet.flow_vol[0].setlb(20.01)
m.fs.strip_sx[1].aqueous_inlet.flow_vol[0].setub(80.01)

# m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H"].setlb(10**-2 * units.g / units.L)
# m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H"].setub(6 * units.g / units.L)

# m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "H"].setub(10 * units.g / units.L)
# m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "H"].setlb(10**-2 * units.g / units.L)

for i in strip_interstage_list:
    m.fs.aq_inter_mixer[i].feed.flow_vol[0].setlb(0.1)
    m.fs.aq_inter_mixer[i].feed.flow_vol[0].setub(3)

# m.fs.load_sx[load_number_of_stages].organic_inlet.conc_mass_comp[0, "DEHPA"].setlb(
#     975.8e3 * 6 / 100
# )
# m.fs.load_sx[load_number_of_stages].organic_inlet.conc_mass_comp[0, "DEHPA"].setub(
#     975.8e3 * 16 / 100
# )


from_json(
    m,
    fname="sx_ascorbic_extraction_circuit_2.json",
    wts=StoreSpec.value(only_not_fixed=True),
)


@m.Objective(sense=maximize)
def objective_function(m):
    # return m.Nd_Ce_recovery
    # return m.Nd_Ce_composition
    return m.tree_recovery


results = solver.solve(m, tee=True)

# print dfs
print("")
print("Degrees of freedom")
print(
    "Ascorbic acid concentration: ",
    m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ascorbic"](),
)
print(
    "Base flowrate: ",
    m.fs.aq_feed_neutral.base_flowrate[0](),
)
for i in load_interstage_list:
    print(
        f"Load interstage {i} flowrate: ",
        m.fs.org_inter_mixer[i].feed.flow_vol[0](),
    )
print(
    "Scrub inlet flowrate: ",
    m.fs.scrub_sx.aqueous_inlet.flow_vol[0](),
)
print(
    "Strip inlet flowrate: ",
    m.fs.strip_sx[1].aqueous_inlet.flow_vol[0](),
)
for i in strip_interstage_list:
    print(
        f"Strip interstage {i} flowrate: ",
        m.fs.aq_inter_mixer[i].feed.flow_vol[0](),
    )


df_tree_recovery_record = pd.DataFrame()
# df_tree_recovery_record = pd.DataFrame(index=vars_of_interest)

df_tree_recovery_record = pd.read_json("tree_recovery.json", orient="split")


def update_tree_dataset(m, df):
    c = df.shape[1] + 1
    df.loc["Loading pH", c] = m.fs.aq_feed_neutral.control_volume.properties_in[
        0.0
    ].pH_phase["liquid"]()
    df.loc["Loading extractant dosage (%v/v)", c] = (
        m.fs.load_sx[load_number_of_stages].organic_inlet.conc_mass_comp[0, "DEHPA"]()
        / 9758
    )
    df.loc["Scrubbing acid (M)", c] = (
        m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H"]() / 1000
    )
    df.loc["Stripping acid (M)", c] = (
        m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "H"]() / 1000
    )
    df.loc["Tree recovery (%)", c] = round(m.tree_recovery() * 100, 3)
    df.loc["Fe concentration (mg/L)", c] = round(
        m.fs.strip_sx[strip_number_of_stages].aqueous_outlet.conc_mass_comp[0, "Fe"](),
        6,
    )
    df.loc["Ascorbic acid concentration (mg/L)", c] = round(
        m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ascorbic"](),
        3,
    )
    df.loc["Base flowrate (L/hr)", c] = round(
        m.fs.aq_feed_neutral.base_flowrate[0](),
        3,
    )
    df.loc["Scrub inlet flowrate (L/hr)", c] = round(
        m.fs.scrub_sx.aqueous_inlet.flow_vol[0](),
        3,
    )
    df.loc["Strip inlet flowrate (L/hr)", c] = round(
        m.fs.strip_sx[1].aqueous_inlet.flow_vol[0](),
        3,
    )
    for i in load_interstage_list:
        df.loc[f"load interstage {i} flowrate (L/hr)", c] = round(
            m.fs.org_inter_mixer[i].feed.flow_vol[0](), 3
        )
    for i in strip_interstage_list:
        df.loc[f"strip interstage {i} flowrate (L/hr)", c] = round(
            m.fs.aq_inter_mixer[i].feed.flow_vol[0](), 3
        )

    return df


# df_tree_recovery_record = update_tree_dataset(m, df_tree_recovery_record)

# df_tree_recovery_record.to_json("tree_recovery.json", orient="split", indent=4)


def update_lb_ub(m, df):

    df.loc["Ascorbic acid concentration (mg/L)", "lb"] = round(
        m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ascorbic"].lb,
        3,
    )
    df.loc["Base flowrate (L/hr)", "lb"] = round(
        m.fs.aq_feed_neutral.base_flowrate[0].lb,
        3,
    )
    df.loc["Scrub inlet flowrate (L/hr)", "lb"] = round(
        m.fs.scrub_sx.aqueous_inlet.flow_vol[0].lb,
        3,
    )
    df.loc["Strip inlet flowrate (L/hr)", "lb"] = round(
        m.fs.strip_sx[1].aqueous_inlet.flow_vol[0].lb,
        3,
    )
    for i in load_interstage_list:
        df.loc[f"load interstage {i} flowrate (L/hr)", "lb"] = round(
            m.fs.org_inter_mixer[i].feed.flow_vol[0].lb, 3
        )
    for i in strip_interstage_list:
        df.loc[f"strip interstage {i} flowrate (L/hr)", "lb"] = round(
            m.fs.aq_inter_mixer[i].feed.flow_vol[0].lb, 3
        )

    df.loc["Ascorbic acid concentration (mg/L)", "ub"] = round(
        m.fs.aq_feed_neutral.inlet.conc_mass_comp[0, "Ascorbic"].ub,
        3,
    )
    df.loc["Base flowrate (L/hr)", "ub"] = round(
        m.fs.aq_feed_neutral.base_flowrate[0](),
        3,
    )
    df.loc["Scrub inlet flowrate (L/hr)", "ub"] = round(
        m.fs.scrub_sx.aqueous_inlet.flow_vol[0].ub,
        3,
    )
    df.loc["Strip inlet flowrate (L/hr)", "ub"] = round(
        m.fs.strip_sx[1].aqueous_inlet.flow_vol[0].ub,
        3,
    )
    for i in load_interstage_list:
        df.loc[f"load interstage {i} flowrate (L/hr)", "ub"] = round(
            m.fs.org_inter_mixer[i].feed.flow_vol[0].ub, 3
        )
    for i in strip_interstage_list:
        df.loc[f"strip interstage {i} flowrate (L/hr)", "ub"] = round(
            m.fs.aq_inter_mixer[i].feed.flow_vol[0].ub, 3
        )

    return df


# row_name = "Tree recovery (%)"
# df_tree_recovery_record_sorted = df_tree_recovery_record.loc[
#     :, df_tree_recovery_record.loc[row_name].sort_values(ascending=False).index
# ]
