#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import (
    ConcreteModel,
    units,
    TransformationFactory,
    Var,
    RangeSet,
    log10,
)
from pyomo.dae.flatten import flatten_dae_components
from pyomo.dae import DerivativeVar
from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)
from idaes.core.util import from_json, StoreSpec

from idaes.core.solvers import get_solver

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution_new import REESolExOgParameters
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package_new_modified import (
    SolventExtractionReactions,
)
import matplotlib.pyplot as plt

# time_duration = 24


def build_model(dosage, number_of_stages, time_duration):
    """
    Method to build a dynamic flowsheet for mixer settler solvent extraction.
    Args:
        dosage: percentage dosage of extractant to the system.
        number_of_stages: number of stages in the mixer settler model.
        time_duration = total time of operation of the model
    Returns:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(
        dynamic=True, time_set=[0, time_duration], time_units=units.hour
    )
    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.reaxn = SolventExtractionReactions()

    m.fs.reaxn.extractant_dosage = dosage

    m.fs.mixer_settler_ex = MixerSettlerExtraction(
        number_of_stages=number_of_stages,
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

    return m


def discretization_scheme(m):
    """
    Discretize the mixer settler solvent extraction model
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None
    """
    m.discretizer = TransformationFactory("dae.collocation")
    m.discretizer.apply_to(m, nfe=12, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")


def copy_first_steady_state(m):
    """
    Function that propagates initial steady state guess to future time points.
    This function is used to initialize all the time discrete variables to the
    initial steady state value.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None
    """
    regular_vars, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
    # Copy initial conditions forward
    for var in time_vars:
        for t in m.fs.time:
            if t == m.fs.time.first():
                continue
            else:
                var[t].value = var[m.fs.time.first()].value


# Fixing inlet conditions


def set_inputs(m, dosage, perturb_time):
    """
    Set inlet conditions to the mixer settler solvent extraction model and fixing
    the parameters of the model.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
        dosage: percentage dosage of extractant to the system.
        perturb_time : time at which a perturbation is added in the flowsheet, should
        be lesser than the time of operation.
    Returns:
        None

    """

    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "H2O"].fix(1e6)
    # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "H"].fix(10.75)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "SO4"].fix(100)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "HSO4"].fix(1e4)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Al"].fix(422.375)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Ca"].fix(109.542)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Cl"].fix(1e-7)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Fe"].fix(688.266)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Sc"].fix(0.032)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Y"].fix(0.124)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "La"].fix(0.986)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Ce"].fix(2.277)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Pr"].fix(0.303)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Nd"].fix(0.946)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Sm"].fix(0.097)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Gd"].fix(0.2584)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Dy"].fix(0.047)

    for t in m.fs.time:
        if t <= perturb_time:
            m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t].fix(62.01)
        else:
            m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t].fix(62.01)
            # m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t].fix(62.01)
        if t <= perturb_time:
            m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Gd"].fix(0.2584)
        else:
            # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Gd"].fix(0.2584)
            m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Gd"].fix(
                0.2584 * 1.5
            )
        if t <= perturb_time * 2:
            m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"].fix(10.75)
        # elif perturb_time <= t < perturb_time * 2:
        #     m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"].fix(8.75)
        else:
            m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"].fix(6.75)
            # m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"].fix(10.75)

    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Kerosene"].fix(820e3)
    # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "DEHPA"].fix(
    #     975.8e3 * dosage / 100
    # )
    for t in m.fs.time:
        if t <= perturb_time * 2:
            m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[t, "DEHPA"].fix(
                975.8e3 * dosage / 100
            )
        else:
            m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[t, "DEHPA"].fix(
                975.8e3 * dosage * 1 / 100
            )
            # m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[t, "DEHPA"].fix(
            #     975.8e3 * dosage / 100
            # )
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Al_o"].fix(1.267e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Ca_o"].fix(2.684e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Fe_o"].fix(2.873e-6)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Sc_o"].fix(1.734)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Y_o"].fix(2.179e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "La_o"].fix(0.000105)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Ce_o"].fix(0.00031)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Pr_o"].fix(3.711e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Nd_o"].fix(0.000165)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Sm_o"].fix(1.701e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Gd_o"].fix(3.357e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Dy_o"].fix(8.008e-6)

    m.fs.mixer_settler_ex.organic_inlet.flow_vol.fix(62.01)

    # Fixing mixer parameters

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
        305.15 * units.K
    )
    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
        305.15 * units.K
    )

    # Fixing settler parameters

    m.fs.mixer_settler_ex.organic_settler[:].unit.area.fix(1)
    m.fs.mixer_settler_ex.aqueous_settler[:].unit.area.fix(1)
    m.fs.mixer_settler_ex.aqueous_settler[:].unit.length.fix(0.1)
    m.fs.mixer_settler_ex.organic_settler[:].unit.length.fix(0.1)


def set_initial_conditions(m):
    """
    Set initial conditions at time=0 for the mixer-settler solvent extraction model
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None

    """

    for e in m.fs.leach_soln.component_list:
        if e not in ["H2O", "HSO4"]:
            m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[
                0, :
            ].conc_mass_comp[e].fix()

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.volume_frac_stream[
        0, :, "aqueous"
    ].fix()
    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[0, :].flow_vol.fix()

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous_inherent_reaction_extent[
        0.0, :, "Ka2"
    ].fix()

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[0, :].flow_vol.fix()
    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[0, :].conc_mass_comp[
        "DEHPA"
    ].fix()

    for e in m.fs.reaxn.element_list:
        m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.heterogeneous_reaction_extent[
            0.0, :, f"{e}_mass_transfer"
        ].fix()

    for e in ["Al", "Ca", "Fe", "Sc"]:
        for s in m.fs.mixer_settler_ex.elements:
            m.fs.mixer_settler_ex.mixer[
                s
            ].unit.mscontactor.heterogeneous_reaction_extent[
                :, :, f"{e}_mass_transfer"
            ].fix(
                0
            )
            m.fs.mixer_settler_ex.mixer[s].unit.mscontactor.organic[
                0.0, 1
            ].conc_mass_comp[f"{e}_o"].fix()
            m.fs.mixer_settler_ex.mixer[s].unit.distribution_extent_constraint[
                :, :, e
            ].deactivate()

    # set variable values in the settler at t=0

    for s in m.fs.mixer_settler_ex.elements:
        for x in m.fs.mixer_settler_ex.aqueous_settler[s].unit.length_domain:
            if x != 0:
                for e in m.fs.leach_soln.component_list:
                    if e not in ["H2O", "HSO4"]:
                        m.fs.mixer_settler_ex.aqueous_settler[s].unit.properties[
                            0, x
                        ].conc_mass_comp[e].fix()
                m.fs.mixer_settler_ex.aqueous_settler[s].unit.properties[
                    0, x
                ].flow_vol.fix()
                m.fs.mixer_settler_ex.aqueous_settler[s].unit.inherent_reaction_extent[
                    0, x, "Ka2"
                ].fix()

        for x in m.fs.mixer_settler_ex.organic_settler[s].unit.length_domain:
            if x != 0:
                for e in m.fs.prop_o.component_list:
                    if e not in ["Kerosene"]:
                        m.fs.mixer_settler_ex.organic_settler[s].unit.properties[
                            0, x
                        ].conc_mass_comp[e].fix()
                m.fs.mixer_settler_ex.organic_settler[s].unit.properties[
                    0, x
                ].flow_vol.fix()


def build_model_and_discretize(dosage, number_of_stages, time_duration):
    """
    Method to build a dynamic model for mixer settler solvent extraction and discretize
    the model.
    Args:
        dosage: percentage dosage of extractant to the system.
        number_of_stages: number of stages in the mixer settler model.
        time_duration = total time of operation of the model
    Returns:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    """

    m = build_model(dosage, number_of_stages, time_duration)

    # m.recovery_Y_integral = Var(m.fs.time)
    # m.differential_Y_recovery = DerivativeVar(m.recovery_Y_integral, wrt=m.fs.time)

    # @m.Constraint(m.fs.time)
    # def integral_recovery(m, t):
    #     return (
    #         m.differential_Y_recovery[t]
    #         == (
    #             1
    #             - (
    #                 m.fs.mixer_settler_ex.aqueous_outlet.conc_mass_comp[t, "Y"]
    #                 * m.fs.mixer_settler_ex.aqueous_outlet.flow_vol[t]
    #             )
    #             / (
    #                 m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Y"]
    #                 * m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t]
    #             )
    #         )
    #         * 100
    #     )

    # m.recovery_Y_integral[0].fix(0)

    # m.actual_recovery_Y = Var(m.fs.time)

    # @m.Constraint(m.fs.time)
    # def actual_recovery(m, t):
    #     if t == 0:
    #         return (
    #             m.actual_recovery_Y[t]
    #             == (
    #                 1
    #                 - (
    #                     m.fs.mixer_settler_ex.aqueous_outlet.conc_mass_comp[0, "Y"]
    #                     * m.fs.mixer_settler_ex.aqueous_outlet.flow_vol[0]
    #                 )
    #                 / (
    #                     m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[0, "Y"]
    #                     * m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[0]
    #                 )
    #             )
    #             * 100
    #         )
    #     else:
    #         return m.actual_recovery_Y[t] * t == m.recovery_Y_integral[t]

    discretization_scheme(m)

    return m


def import_steady_value(m, path_name):
    """
    A function to import the steady state values of the mixer-settler solvent extraction
    model to the dynamic model for initializing it.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
        path_name: name of the path of the json file
    Returns:
        None
    """
    from_json(m, fname=path_name, wts=StoreSpec.value())


def initialize_set_input_and_initial_conditions(m, dosage, perturb_time):
    """
    Function to initialize, set inlet values and give initial conditions to the dynamic
    mixer settler solvent extraction model.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
        dosage: percentage dosage of extractant to the system.
        perturb_time : time at which a perturbation is added in the flowsheet, should
        be lesser than the time of operation.
    Returns:
        None

    """

    copy_first_steady_state(m)
    set_inputs(m, dosage, perturb_time)
    set_initial_conditions(m)


def solve_model(m):
    """
    A function to solve the initialized mixer-settler solvent extraction model.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None
    """
    solver = get_solver("ipopt_v2")
    results = solver.solve(m, tee=True)
    return results


def main(dosage, number_of_stages, time_duration, perturb_time, path_name):
    """
    Function to build a dynamic model, discretize it, initialize it, set input values
    and give initial conditions, then solve the model.
    Args:
        dosage: percentage dosage of extractant to the system.
        number_of_stages: number of stages in the mixer settler model.
        time_duration = total time of operation of the model
        path_name: name of the path of the json file
        perturb_time : time at which a perturbation is added in the flowsheet, should
        be lesser than the time of operation.
    Returns:
        m: ConcreteModel object with the mixer-settler solvent extraction system.

    """
    m = build_model_and_discretize(dosage, number_of_stages, time_duration)
    import_steady_value(m, path_name)
    initialize_set_input_and_initial_conditions(m, dosage, perturb_time)
    results = solve_model(m)

    return m, results


dosage = 5
number_of_stages = 3
time_duration = 24
perturb_time = 4
if __name__ == "__main__":
    m, results = main(
        dosage,
        number_of_stages,
        time_duration,
        perturb_time,
        path_name="mixer_settler_extraction.json",
    )

percentage_recovery = {}

# for e in m.fs.leach_soln.component_list:
#     if e not in ["H2O", "H", "SO4", "HSO4", "Cl"]:
#         percentage_recovery[e] = [
#             (
#                 1
#                 - (
#                     m.fs.mixer_settler_ex.aqueous_outlet.conc_mass_comp[t, e]()
#                     * m.fs.mixer_settler_ex.aqueous_outlet.flow_vol[t]()
#                 )
#                 / (
#                     m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, e]()
#                     * m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t]()
#                 )
#             )
#             * 100
#             for t in m.fs.time
#         ]

for e in m.fs.leach_soln.component_list:
    if e not in ["H2O", "H", "SO4", "HSO4", "Cl"]:
        percentage_recovery[e] = [
            (
                (
                    m.fs.mixer_settler_ex.organic_outlet.conc_mass_comp[t, f"{e}_o"]()
                    * m.fs.mixer_settler_ex.organic_outlet.flow_vol[t]()
                    - m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[t, f"{e}_o"]()
                    * m.fs.mixer_settler_ex.organic_inlet.flow_vol[t]()
                )
                / (
                    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, e]()
                    * m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t]()
                )
            )
            * 100
            for t in m.fs.time
        ]


Gd_recovery_fraction = [
    (
        m.fs.mixer_settler_ex.organic_outlet.conc_mass_comp[t, "Gd_o"]()
        * m.fs.mixer_settler_ex.organic_outlet.flow_vol[t]()
    )
    / (
        m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Gd"]()
        * m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t]()
    )
    for t in m.fs.time
]


Gd_loss_fraction = [
    (
        m.fs.mixer_settler_ex.aqueous_outlet.conc_mass_comp[t, "Gd"]()
        * m.fs.mixer_settler_ex.aqueous_outlet.flow_vol[t]()
    )
    / (
        m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Gd"]()
        * m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t]()
    )
    for t in m.fs.time
]

Gd_accum_fraction = [
    sum(
        m.fs.mixer_settler_ex.mixer[s].unit.mscontactor.aqueous_material_accumulation[
            t, 1, "liquid", "Gd"
        ]()
        + m.fs.mixer_settler_ex.mixer[s].unit.mscontactor.organic_material_accumulation[
            t, 1, "organic", "Gd_o"
        ]()
        for s in m.fs.mixer_settler_ex.elements
    )
    / (
        m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "Gd"]()
        * m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t]()
        / (m.fs.leach_soln.mw["Gd"]() * 1e6)
    )
    for t in m.fs.time
]

plt.rcParams.update(
    {
        "figure.max_open_warning": 0,
        "figure.dpi": 300,
        "figure.titlesize": 16,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "axes.linewidth": 2,
        "lines.linewidth": 2,
        "lines.markersize": 8,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "savefig.bbox": "tight",
        "legend.fontsize": "large",
    }
)

# REE_list = []
# for e in m.fs.leach_soln.component_list:
#     if e in ["Y", "Dy", "Gd", "La"]:
#         REE_list.append(e)
#         plt.plot(
#             m.fs.time,
#             percentage_recovery[e],
#         )
# plt.legend(REE_list)

# plt.show()

# for s in RangeSet(number_of_stages):
#     plt.plot(
#         m.fs.time,
#         m.fs.mixer_settler_ex.mixer[s]
#         .unit.mscontactor.organic[:, 1]
#         .conc_mass_comp["Y_o"](),
#     )
# plt.legend(["stage 1", "stage 2", "stage 3"])


# fig, ax = plt.subplots(1, 3, figsize=(7, 9), dpi=300)

# fig.suptitle("pH perturbation effect on Gd")
# ax[0].step(
#     m.fs.time,
#     [
#         -log10(m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"]() / 1000)
#         for t in m.fs.time
#     ],
#     linewidth=3,
# )
# ax[0].set_xlabel("Time, hr")
# ax[0].set_ylabel("pH")
# ax[0].set_title("Aqueous feed pH")
# ax[0].axvline(4, linestyle="--", color="green", linewidth=2)
# ax[1].plot(
#     m.fs.time,
#     m.fs.mixer_settler_ex.organic_settler[1]
#     .unit.properties[:, 1]
#     .conc_mass_comp["Gd_o"](),
#     linewidth=3,
#     label="stage 1",
# )
# ax[1].set_xlabel("Time, hr")
# ax[1].set_ylabel("Stage 1, Concentration, mg/L")
# ax2 = ax[1].twinx()
# ax2.plot(
#     m.fs.time,
#     m.fs.mixer_settler_ex.organic_settler[3]
#     .unit.properties[:, 1]
#     .conc_mass_comp["Gd_o"](),
#     linewidth=3,
#     color="red",
#     label="stage 3",
# )
# ax[1].axvline(4, linestyle="--", color="green", linewidth=2)
# ax2.set_ylabel("Stage 3, Concentration, mg/L")
# ax[1].set_title("Stage 1 and 3 organic settler outlet concentration profile")
# ax[1].ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
# ax2.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
# handles1, labels1 = ax[1].get_legend_handles_labels()
# handles2, labels2 = ax2.get_legend_handles_labels()
# all_handles = handles1 + handles2
# all_labels = labels1 + labels2
# ax[1].legend(all_handles, all_labels, loc="upper left")
# ax[2].plot(m.fs.time, percentage_recovery["Gd"], linewidth=3)
# ax[2].axvline(4, linestyle="--", color="green", linewidth=2)
# ax[2].set_xlabel("Time, hr")
# ax[2].set_ylabel("Recovery %")
# ax[2].set_title("Percentage recovery profile")
# plt.tight_layout()

# # Add centered labels below each subplot (small font)
# plt.subplots_adjust(bottom=0.10)
# ax[0].text(
#     0.5,
#     -0.25,
#     "(a) Aqueous feed pH",
#     transform=ax[0].transAxes,
#     ha="center",
#     fontsize=8,
# )
# ax[1].text(
#     0.5,
#     -0.25,
#     "(b) Organic settler outlet concentration (stages 1 & 3)",
#     transform=ax[1].transAxes,
#     ha="center",
#     fontsize=8,
# )
# ax[2].text(
#     0.5,
#     -0.25,
#     "(c) Percentage recovery",
#     transform=ax[2].transAxes,
#     ha="center",
#     fontsize=8,
# )

# plt.show()

fig, ax = plt.subplots(1, 3, figsize=(15, 4), dpi=300)

ax[0].step(
    m.fs.time,
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Gd"](),
    linewidth=3,
    color="black",
    label="Concentration",
)
ax[0].axvline(perturb_time, linestyle="--", color="black", linewidth=2)
ax[0].set_xlabel("Time, hr")
ax[0].set_ylabel("Concentration, mg/L")
ax[0].set_ylim(0.225, 0.425)
# ax[0].step(
#     m.fs.time,
#     m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[:](),
#     linewidth=3,
#     color="black",
#     label="Flowrate",
# )
# ax[0].axvline(perturb_time, linestyle="--", color="black", linewidth=2)
# ax[0].set_xlabel("Time, hr")
# ax[0].set_ylabel("Flowrate L/hr")
# ax[0].set_ylim(60, 70)

ax_0 = ax[0].twinx()

ax_0.step(
    m.fs.time,
    [
        -log10(m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"]() / 1000)
        for t in m.fs.time
    ],
    linewidth=3,
    color="brown",
    label="pH",
)
# ax_0.step(
#     m.fs.time,
#     [
#         m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[t, "DEHPA"]() * 100 / 975.8e3
#         for t in m.fs.time
#     ],
#     linewidth=3,
#     label="Extractant dosage",
#     color="red",
# )
ax_0.set_xlabel("Time, hr")
ax_0.set_ylabel("pH")
# ax_0.set_ylabel("Dosage % v/v")
ax_0.set_ylim(1.9, 2.2)
ax_0.axvline(perturb_time * 2, linestyle="--", color="brown", linewidth=2)
handles0, labels0 = ax[0].get_legend_handles_labels()
handles0_twin, labels0_twin = ax_0.get_legend_handles_labels()
ax[0].legend(handles0 + handles0_twin, labels0 + labels0_twin, loc="lower right")


ax_1 = ax[1].twinx()
ax[1].plot(m.fs.time, Gd_recovery_fraction, linewidth=3, label="Recovery")
ax_1.plot(m.fs.time, Gd_loss_fraction, linewidth=3, label="Loss", color="purple")
ax[1].axvline(perturb_time, linestyle="--", color="black", linewidth=2)
ax[1].axvline(perturb_time * 2, linestyle="--", color="brown", linewidth=2)
ax[1].set_xlabel("Time, hr")
ax[1].set_ylabel("Recovery fraction")
ax[1].set_ylim(0.22, 0.37)
ax_1.set_ylim(0.62, 0.67)
ax_1.set_ylabel("Loss fraction")
handles1, labels1 = ax[1].get_legend_handles_labels()
handles1_twin, labels1_twin = ax_1.get_legend_handles_labels()
ax[1].legend(handles1 + handles1_twin, labels1 + labels1_twin, loc="lower right")

ax[2].plot(m.fs.time, Gd_accum_fraction, linewidth=3, label="Recovery", color="green")
ax[2].axvline(perturb_time, linestyle="--", color="black", linewidth=2)
ax[2].axvline(perturb_time * 2, linestyle="--", color="brown", linewidth=2)
ax[2].set_xlabel("Time, hr")
ax[2].set_ylabel("Accumulation fraction")


ax[0].set_axisbelow(True)  # Forces gridlines behind bars/plots
ax[0].grid(True)
ax[1].set_axisbelow(True)  # Forces gridlines behind bars/plots
ax[1].grid(True)
ax[2].set_axisbelow(True)  # Forces gridlines behind bars/plots
ax[2].grid(True)


plt.tight_layout()


# fig.suptitle("Aqueous feed flowrate perturbation effect on Gd")
# ax[0].plot(
#     m.fs.time,
#     m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[:](),
#     linewidth=3,
# )
# ax[0].axvline(4, linestyle="--", color="green", linewidth=2)
# ax[0].set_xlabel("Time, hr")
# ax[0].set_ylabel("Flowrate L/hr")
# ax[0].set_title("Aqueous feed flowrate")
# ax[1].plot(
#     m.fs.time,
#     m.fs.mixer_settler_ex.organic_settler[1]
#     .unit.properties[:, 1]
#     .conc_mass_comp["Gd_o"](),
#     linewidth=3,
#     label="stage 1",
# )
# ax[1].axvline(4, linestyle="--", color="green", linewidth=2)
# ax[1].set_xlabel("Time, hr")
# ax[1].set_ylabel("Stage 1, Concentration, mg/L")
# ax2 = ax[1].twinx()
# ax2.plot(
#     m.fs.time,
#     m.fs.mixer_settler_ex.organic_settler[3]
#     .unit.properties[:, 1]
#     .conc_mass_comp["Gd_o"](),
#     linewidth=3,
#     color="red",
#     label="stage 3",
# )
# ax2.set_ylabel("Stage 3, Concentration, mg/L")
# ax[1].set_title("Stage 1 and 3 organic settler outlet concentration profile")
# handles1, labels1 = ax[1].get_legend_handles_labels()
# handles2, labels2 = ax2.get_legend_handles_labels()
# all_handles = handles1 + handles2
# all_labels = labels1 + labels2
# ax[1].legend(all_handles, all_labels, loc="upper left")
# ax[1].ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
# ax2.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
# ax[2].plot(m.fs.time, percentage_recovery["Gd"], linewidth=3)
# ax[2].axvline(4, linestyle="--", color="green", linewidth=2)
# ax[2].set_xlabel("Time, hr")
# ax[2].set_ylabel("Recovery %")
# ax[2].set_title("Percentage recovery profile")
# plt.tight_layout()
# plt.subplots_adjust(bottom=0.10)
# ax[0].text(
#     0.5,
#     -0.25,
#     "(d) Aqueous feed flowrate",
#     transform=ax[0].transAxes,
#     ha="center",
#     fontsize=8,
# )
# ax[1].text(
#     0.5,
#     -0.25,
#     "(e) Organic settler outlet concentration (stages 1 & 3)",
#     transform=ax[1].transAxes,
#     ha="center",
#     fontsize=8,
# )
# ax[2].text(
#     0.5,
#     -0.25,
#     "(f) Percentage recovery",
#     transform=ax[2].transAxes,
#     ha="center",
#     fontsize=8,
# )


# fig, ax = plt.subplots(1, 3, figsize=(18, 5), dpi=300)

# fig.suptitle("pH perturbation effect on Gd")
# ax[0].step(
#     m.fs.time,
#     [
#         -log10(m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[t, "H"]() / 1000)
#         for t in m.fs.time
#     ],
#     linewidth=3,
# )
# ax[0].set_xlabel("Time, hr")
# ax[0].set_ylabel("pH")
# ax[0].set_ylim(1.9, 2.3)
# ax[0].set_title("Aqueous feed pH")
# ax[0].axvline(4, linestyle="--", color="green", linewidth=2)
# ax[1].plot(
#     m.fs.time,
#     m.fs.mixer_settler_ex.organic_settler[1]
#     .unit.properties[:, 1]
#     .conc_mass_comp["Gd_o"](),
#     linewidth=3,
#     label="stage 1",
# )
# ax[1].set_xlabel("Time, hr")
# ax[1].set_ylabel("Stage 1, Conc. (mg/L)")
# ax2 = ax[1].twinx()
# ax2.plot(
#     m.fs.time,
#     m.fs.mixer_settler_ex.organic_settler[3]
#     .unit.properties[:, 1]
#     .conc_mass_comp["Gd_o"](),
#     linewidth=3,
#     color="red",
#     label="stage 3",
# )
# ax[1].axvline(4, linestyle="--", color="green", linewidth=2)
# ax2.set_ylabel("Stage 3, Conc. (mg/L)")
# ax[1].set_title("Organic settler outlet concentration profile")
# ax[1].ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
# ax2.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
# handles1, labels1 = ax[1].get_legend_handles_labels()
# handles2, labels2 = ax2.get_legend_handles_labels()
# all_handles = handles1 + handles2
# all_labels = labels1 + labels2
# ax[1].legend(all_handles, all_labels, loc="lower right")
# ax[2].plot(m.fs.time, percentage_recovery["Gd"], linewidth=3)
# ax[2].axvline(4, linestyle="--", color="green", linewidth=2)
# ax[2].set_xlabel("Time, hr")
# ax[2].set_ylabel("Recovery %")
# ax[2].set_title("Percentage recovery profile")
# ax[0].text(
#     4,
#     2,
#     " perturbation",
#     fontsize=12,
#     va="top",
#     ha="left",
#     color="black",
# )
# ax[2].text(
#     4,
#     26,
#     " perturbation",
#     fontsize=12,
#     va="top",
#     ha="left",
#     color="black",
# )
# ax[1].text(
#     4,
#     0.067,
#     " perturbation",
#     fontsize=12,
#     va="top",
#     ha="left",
#     color="black",
# )
# ax[2].set_ylim(24, 27.5)

# fig.subplots_adjust(wspace=1.5, bottom=0.12)

# # Add centered labels below each subplot (small font)
# plt.tight_layout()
# plt.show()
