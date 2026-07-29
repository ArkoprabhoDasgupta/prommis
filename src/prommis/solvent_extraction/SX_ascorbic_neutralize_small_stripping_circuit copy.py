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
)


m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()


# define stages
dosage = 20
strip_number_of_stages = 2
strip_stage_list = RangeSet(1, strip_number_of_stages)
strip_interstage_list = RangeSet(1, strip_number_of_stages - 1)

m.fs.reaxn.extractant_dosage = dosage

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


m.organic_scrub_to_strip = Arc(
    source=m.fs.scrub_sx.organic_outlet,
    destination=m.fs.strip_sx[strip_number_of_stages].organic_inlet,
)


TransformationFactory("network.expand_arcs").apply_to(m)


# define inlet conditions


# define aqueous interstage inlet
for i in strip_interstage_list:
    for e in m.fs.leach_soln.component_list:
        if e not in ["H2O", "H", "Cl"]:
            m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, e].fix(1e-7)
    m.fs.aq_inter_mixer[i].feed.flow_vol.fix(1)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H2O"].fix(1e6)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "H"].fix(1 * units.g / units.L)
    m.fs.aq_inter_mixer[i].feed.conc_mass_comp[0, "Cl"].fix(
        1 * 35.5 * units.g / units.L
    )

m.fs.aq_inter_mixer[:].mixed_state[0.0].pressure.fix(101235)
m.fs.aq_inter_mixer[:].mixed_state[0.0].temperature.fix(303.5)

# define scrub sx aqueous inlet
m.fs.scrub_sx.aqueous_inlet.flow_vol.fix(20.01)
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "Cl"]:
        m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, e].fix(1e-7)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "H"].fix(0.1 * units.g / units.L)
m.fs.scrub_sx.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(0.1 * 35.5 * units.g / units.L)

m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Ce_o"].fix(0.635)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "DEHPA"].fix(77952.034)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Dy_o"].fix(0.083)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(5.573)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820000)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "La_o"].fix(0.2277)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Nd_o"].fix(0.2705)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Pr_o"].fix(0.1607)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Sm_o"].fix(0.059)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Y_o"].fix(0.319)
m.fs.scrub_sx.organic_inlet.conc_mass_comp[0, "Gd_o"].fix(0.3559)
m.fs.scrub_sx.organic_inlet.flow_vol.fix(65.01)

m.fs.scrub_sx.mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.scrub_sx.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
m.fs.scrub_sx.mscontactor.organic[:, :].temperature.fix(305.15 * units.K)


# define strip sx aqueous inlet
m.fs.strip_sx[1].aqueous_inlet.flow_vol.fix(80.01)
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "Cl"]:
        m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, e].fix(1e-7)
m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "H"].fix(12 * units.g / units.L)
m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[0, "Cl"].fix(6 * 35.5 * units.g / units.L)
m.fs.strip_sx[:].mscontactor.volume[:].fix(0.4 * units.m**3)
m.fs.strip_sx[:].mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
m.fs.strip_sx[:].mscontactor.organic[:, :].temperature.fix(305.15 * units.K)


solver = get_solver("ipopt_v2")
solver.options["halt_on_ampl_error"] = "yes"
solver.options["max_iter"] = 1000

print(dof(m))


sx_tear_guesses = {
    "flow_vol": {0: 60.01},
    "temperature": {0: 303.15},
    "pressure": {0: 101325.0},
    "conc_mass_comp": {
        (0, "H2O"): 1e6,
        (0, "H"): 6e3,
        (0, "Cl"): 35.5 * 6e3,
        (0, "Ascorbic"): 1e-9,
        (0, "SO4"): 1e-9,
        (0, "HSO4"): 1e-9,
        (0, "La"): 1e-9,
        (0, "Ce"): 1e-9,
        (0, "Pr"): 1e-9,
        (0, "Nd"): 1e-9,
        (0, "Sm"): 1e-9,
        (0, "Gd"): 1e-9,
        (0, "Dy"): 1e-9,
        (0, "Y"): 1e-9,
        (0, "Fe"): 1e-9,
    },
}


sx_initializer = SolventExtractionInitializer()
scrub_sx_units = m.fs.scrub_sx
strip_sx_units = [m.fs.strip_sx[s] for s in strip_stage_list]

mixer_initializer = MixerInitializer()
strip_interstage_mixer = [m.fs.aq_inter_mixer[s] for s in strip_interstage_list]
initialized_mixers = set()


sx_initializer.initialize(scrub_sx_units)
# assert 1 == 2

# # 5. Initialize the scrub and strip units before the decomposition loop
# for unit in [m.fs.scrub_sx, *strip_sx_units]:
#     print(f"Initializing {unit}")
#     sx_initializer.initialize(unit)


def function(unit):
    if unit in strip_interstage_mixer and unit not in initialized_mixers:
        print(f"Initializing mixer {unit}")
        mixer_initializer.initialize(unit)
        initialized_mixers.add(unit)
    if unit in [m.fs.scrub_sx, *strip_sx_units]:
        print(f"Initializing {unit}")
        sx_initializer.initialize(unit)


seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 3

# Using the SD tool
G = seq.create_graph(m)
heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
order = seq.calculation_order(G)

for o in heuristic_tear_set:
    print(o.name)

# for o in order:
#     print(o[0].name)

# assert 1 == 2

seq.set_guesses_for(m.fs.aq_inter_mixer[1].sx, sx_tear_guesses)
# seq.set_guesses_for(m.fs.aq_inter_mixer[2].sx, sx_tear_guesses)
# seq.set_guesses_for(m.fs.aq_inter_mixer[3].sx, sx_tear_guesses)

seq.run(m, function)

results = solver.solve(m, tee=True)

percentage_recovery = {}
for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic"]:
        percentage_recovery[e] = [
            (
                (
                    (
                        m.fs.strip_sx[
                            strip_number_of_stages
                        ].aqueous_outlet.conc_mass_comp[t, e]()
                        * m.fs.strip_sx[strip_number_of_stages].aqueous_outlet.flow_vol[
                            t
                        ]()
                        - m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[t, e]()
                        * m.fs.strip_sx[1].aqueous_inlet.flow_vol[t]()
                    )
                )
                / (
                    m.fs.scrub_sx.organic_inlet.conc_mass_comp[t, f"{e}_o"]()
                    * m.fs.scrub_sx.organic_inlet.flow_vol[t]()
                )
            )
            * 100
            for t in m.fs.time
        ]


percentage_recovery["tree"] = [
    (
        (
            (
                sum(
                    m.fs.strip_sx[strip_number_of_stages].aqueous_outlet.conc_mass_comp[
                        t, e
                    ]()
                    * m.fs.strip_sx[strip_number_of_stages].aqueous_outlet.flow_vol[t]()
                    for e in m.fs.leach_soln.component_list
                    if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
                )
                - sum(
                    m.fs.strip_sx[1].aqueous_inlet.conc_mass_comp[t, e]()
                    * m.fs.strip_sx[1].aqueous_inlet.flow_vol[t]()
                    for e in m.fs.leach_soln.component_list
                    if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
                )
            )
        )
        / sum(
            m.fs.scrub_sx.organic_inlet.conc_mass_comp[t, f"{e}_o"]()
            * m.fs.scrub_sx.organic_inlet.flow_vol[t]()
            for e in m.fs.leach_soln.component_list
            if e not in ["H2O", "H", "SO4", "HSO4", "Cl", "Ascorbic", "Fe"]
        )
    )
    * 100
    for t in m.fs.time
]


# print(percentage_recovery)
