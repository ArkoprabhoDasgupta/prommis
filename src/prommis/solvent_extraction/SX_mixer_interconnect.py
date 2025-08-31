from pyomo.environ import (
    ConcreteModel,
    units,
    RangeSet,
    TransformationFactory,
    Var,
    maximize,
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
from prommis.solvent_extraction.solvent_extraction_reaction_package_new import (
    SolventExtractionReactions,
)

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

dosage = 5

number_of_stages = 2
number_of_interstage_mixers = number_of_stages - 1

stage_list = RangeSet(1, number_of_stages)
interstage_list = RangeSet(1, number_of_interstage_mixers)

# m.fs.reaxn.extractant_dosage = dosage

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

m.fs.interstage_mixer = Mixer(
    interstage_list,
    property_package=m.fs.prop_o,
    num_inlets=2,
    inlet_list=["organic_out", "feed"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

for i in stage_list:
    if i != number_of_stages:
        m.add_component(
            f"aqueous_sx_{i}_to_{i+1}",
            Arc(
                source=m.fs.mixer_settler_sx[i].aqueous_outlet,
                destination=m.fs.mixer_settler_sx[i + 1].aqueous_inlet,
            ),
        )
        m.add_component(
            f"interstage_mixer_{i}_to_organic_sx_{i}",
            Arc(
                source=m.fs.interstage_mixer[i].outlet,
                destination=m.fs.mixer_settler_sx[i].organic_inlet,
            ),
        )
    if i != 1:
        m.add_component(
            f"organic_sx_{i}_to_interstage_mixer_{i-1}",
            Arc(
                source=m.fs.mixer_settler_sx[i].organic_outlet,
                destination=m.fs.interstage_mixer[i - 1].organic_out,
            ),
        )

TransformationFactory("network.expand_arcs").apply_to(m)

m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "H"].fix(10.75)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "SO4"].fix(100)
m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e4)
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

m.fs.mixer_settler_sx[number_of_stages].organic_inlet.extractant_dosage.fix(5)
m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "Kerosene"].fix(
    820e3
)
# m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
#     975.8e3 * dosage / 100
# )
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

m.fs.mixer_settler_sx[:].mixer[:].unit.mscontactor.volume[:].fix(1e-4 * units.m**3)

m.fs.mixer_settler_sx[:].organic_settler[:].unit.area.fix(1e-2)
m.fs.mixer_settler_sx[:].aqueous_settler[:].unit.area.fix(1e-2)
m.fs.mixer_settler_sx[:].aqueous_settler[:].unit.length.fix(1e-2)
m.fs.mixer_settler_sx[:].organic_settler[:].unit.length.fix(1e-2)

m.fs.interstage_mixer[1].feed.extractant_dosage.fix(3)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Kerosene"].fix(820e3)
# m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "DEHPA"].fix(975.8e3 * dosage / 100)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Al_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Ca_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Fe_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Sc_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Y_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "La_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Ce_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Pr_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Nd_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Sm_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Gd_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.conc_mass_comp[0, "Dy_o"].fix(1e-9)
m.fs.interstage_mixer[1].feed.flow_vol.fix(12.01)

m.fs.interstage_mixer[1].mixed_state[0.0].temperature.fix(305.15 * units.K)
m.fs.interstage_mixer[1].mixed_state[0.0].pressure.fix(1e5 * units.Pa)

m.fs.mixer_settler_sx[:].mixer_organic_settler_arc_1_expanded.conc_mass_comp_equality[
    0, "DEHPA"
].deactivate()

m.interstage_mixer_1_to_organic_sx_1_expanded.conc_mass_comp_equality[
    0, "DEHPA"
].deactivate()
m.organic_sx_2_to_interstage_mixer_1_expanded.conc_mass_comp_equality[
    0, "DEHPA"
].deactivate()

trial_element_list = ["Y", "Dy", "Gd", "Sm", "Nd", "Ce"]

m.percentage_recovery = Var(trial_element_list, initialize=1, bounds=(0, 100))


@m.Constraint(trial_element_list)
def recovery_constraint(m, e):
    return (
        m.percentage_recovery[e]
        == (
            (
                m.fs.mixer_settler_sx[1].organic_outlet.conc_mass_comp[0, f"{e}_o"]
                * m.fs.mixer_settler_sx[1].organic_outlet.flow_vol[0]
                - m.fs.mixer_settler_sx[number_of_stages].organic_inlet.conc_mass_comp[
                    0, f"{e}_o"
                ]
                * m.fs.mixer_settler_sx[number_of_stages].organic_inlet.flow_vol[0]
            )
            / (
                m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, f"{e}"]
                * m.fs.mixer_settler_sx[1].aqueous_inlet.flow_vol[0]
            )
        )
        * 100
    )


@m.Constraint()
def inlet_constraint(m):
    return (
        m.fs.mixer_settler_sx[1].aqueous_inlet.conc_mass_comp[0, "H"]
        >= 0.1 * units.mg / units.L
    )


# @m.Constraint()
# def Dy_constraint(m):
#     return m.percentage_recovery["Dy"] <= 35


# @m.Constraint()
# def Gd_constraint(m):
#     return m.percentage_recovery["Gd"] <= 50


# @m.Objective(sense=maximize)
# def objective_function(m):
#     return m.percentage_recovery["Gd"]


print(degrees_of_freedom(m))

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 100000
solver.solve(m, tee=True)
