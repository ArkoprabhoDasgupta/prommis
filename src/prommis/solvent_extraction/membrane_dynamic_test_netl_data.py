from pyomo.environ import (
    ConcreteModel,
    units,
    Suffix,
    TransformationFactory,
    log10,
    Var,
)
from pyomo.dae.flatten import flatten_dae_components

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.util.model_statistics import degrees_of_freedom as dof

from prommis.solvent_extraction.membrane_module_property_package import (
    MembraneSXModuleParameters,
)

from prommis.solvent_extraction.membrane_channel_property_package import (
    MembraneSXChannelParameters,
)
from idaes.core.util import from_json


from prommis.solvent_extraction.membrane_solvent_extraction import (
    MembraneSolventExtraction,
    MembraneSolventExtractionInitializer,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from idaes.core.util import DiagnosticsToolbox

m = ConcreteModel()

time_duration = 12

m.fs = FlowsheetBlock(dynamic=True, time_set=[0, time_duration], time_units=units.hour)

# m.fs = FlowsheetBlock(dynamic=True)

m.fs.mem_prop = MembraneSXModuleParameters()
m.fs.mem_channel = MembraneSXChannelParameters()
# m.fs.leach_soln = LeachSolutionParameters()
m.fs.mem_prop.extractant_dosage = 10

m.fs.membrane_module = MembraneSolventExtraction(
    feed_phase={
        "property_package": m.fs.mem_channel,
        "flow_direction": FlowDirection.forward,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "energy_balance_type": EnergyBalanceType.none,
    },
    strip_phase={
        "property_package": m.fs.mem_channel,
        "flow_direction": FlowDirection.backward,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "energy_balance_type": EnergyBalanceType.none,
    },
    membrane_phase={
        "property_package": m.fs.mem_prop,
    },
    finite_elements=5,
    transformation_method="dae.finite_difference",
    transformation_scheme="BACKWARD",
    # collocation_points=2,
    tube_inner_radius=100 * units.micrometer,
    tube_outer_radius=150 * units.micrometer,
    shell_radius=(67 / 2) * units.mm,
    number_of_tubes=6300,
)

m.fs.membrane_module.module_length.fix(0.254 * 2 * units.m)

m.discretizer = TransformationFactory("dae.finite_difference")
m.discretizer.apply_to(m, nfe=5, wrt=m.fs.time, scheme="BACKWARD")


# def copy_first_steady_state(m):
#     """
#     Function that propagates initial steady state guess to future time points.
#     This function is used to initialize all the time discrete variables to the
#     initial steady state value.
#     """
#     regular_vars, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
#     # Copy initial conditions forward
#     for var in time_vars:
#         for t in m.fs.time:
#             if t == m.fs.time.first():
#                 continue
#             else:
#                 var[t].value = var[m.fs.time.first()].value


# from_json(m, fname="membrane_solvent_extraction.json")
# copy_first_steady_state(m)

# Feed phase inlet conditions

m.fs.membrane_module.feed_phase_inlet.flow_vol.fix(17 * units.mL / units.min)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Ca"].fix(
    1535976 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Al"].fix(
    310800 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Fe"].fix(
    131300 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "La"].fix(
    273 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Ce"].fix(
    573 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Pr"].fix(
    70 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Nd"].fix(
    315 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Sm"].fix(
    65.9 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Gd"].fix(
    67.5 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Dy"].fix(
    52.4 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Y"].fix(
    231 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Sc"].fix(
    20 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "H"].fix(
    10**-2.06 * 1 * units.gram / units.L
)
# m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "H"].fix(
#     0.01 * units.gram / units.L
# )
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "H2O"].fix(1e6)

# Strip phase inlet conditions

m.fs.membrane_module.strip_phase_inlet.flow_vol.fix(54 * units.mL / units.min)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Ca"].fix(
    4081 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Al"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Fe"].fix(
    515 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "La"].fix(
    1.89 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Ce"].fix(
    3.19 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Pr"].fix(
    0.398 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Nd"].fix(
    1.24 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Sm"].fix(
    0.15 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Gd"].fix(
    0.157 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Dy"].fix(
    0.119 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Y"].fix(
    0.6 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Sc"].fix(
    0.649 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "H"].fix(
    3 * 1 * units.gram / units.L
)
# m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "H"].fix(
#     0.5 * units.gram / units.L
# )
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "H2O"].fix(1e6)

for e in m.fs.mem_channel.component_list:
    if e not in ["H2O"]:
        m.fs.membrane_module.feed_phase.properties[0.0, :].conc_mass_comp[e].fix(1e-10)
        m.fs.membrane_module.feed_phase.properties[0.0, :].flow_vol.fix(17)

        m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp[e].fix(1e-10)
        m.fs.membrane_module.strip_phase.properties[0.0, :].flow_vol.fix(54)

# for e in m.fs.mem_prop.component_list:
#     m.fs.membrane_module.conc_mol_membrane_comp[0, 0.0, :, e].fix()

print("Degrees of freedom:", dof(m))

m.scaling_factor = Suffix(direction=Suffix.EXPORT)

for t in m.fs.time:
    for z in m.fs.membrane_module.feed_phase.length_domain:
        set_scaling_factor(
            m.fs.membrane_module.feed_phase.properties[t, z].conc_mass_comp["H"],
            1,
        )
        set_scaling_factor(
            m.fs.membrane_module.feed_phase.properties[t, z].pH_constraint,
            1,
        )
        for e in m.fs.mem_prop.component_list:
            set_scaling_factor(
                m.fs.membrane_module.feed_phase.properties[t, z].conc_mol_comp[e],
                1e4,
            )
            set_scaling_factor(
                m.fs.membrane_module.feed_phase.properties[t, z].conc_mass_comp[e],
                1,
            )
            set_scaling_factor(
                m.fs.membrane_module.strip_phase.properties[t, z].conc_mol_comp[e],
                1e4,
            )
            set_scaling_factor(
                m.fs.membrane_module.strip_phase.properties[t, z].conc_mass_comp[e],
                1,
            )

            for r in m.fs.membrane_module.r:
                set_scaling_factor(
                    m.fs.membrane_module.conc_mol_membrane_comp[t, z, r, e],
                    1e4,
                )

scaling = TransformationFactory("core.scale_model")
scaled_model = scaling.create_using(m, rename=False)

# # initializer = MembraneSolventExtractionInitializer()
# # initializer.initialize(scaled_model.fs.membrane_module)
# # # initializer.initialize(m.fs.membrane_module)

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 10000
# results = solver.solve(m, tee=True)
results = solver.solve(scaled_model, tee=True)

# # dt = DiagnosticsToolbox(scaled_model)

# scaling.propagate_solution(scaled_model, m)

# # import pandas as pd

# # index_list = []
# # for e in m.fs.mem_channel.component_list:
# #     if e == "H2O":
# #         pass
# #     elif e == "H":
# #         index_list.append(f"{e}_feed")
# #         index_list.append("pH_feed")
# #         index_list.append("model_pH_feed")
# #         index_list.append(f"{e}_strip")
# #         index_list.append("pH_strip")
# #         index_list.append("model_pH_strip")
# #     else:
# #         index_list.append(f"{e}_feed")
# #         index_list.append(f"{e}_strip")

# # df = pd.DataFrame(
# #     index=index_list,
# #     columns=m.fs.membrane_module.feed_phase.length_domain,
# # )

# # for z in m.fs.membrane_module.feed_phase.length_domain:
# #     for e in m.fs.mem_channel.component_list:
# #         if e == "H":
# #             df.loc[f"{e}_feed", z] = m.fs.membrane_module.feed_phase.properties[
# #                 0, z
# #             ].conc_mol_comp[e]()
# #             df.loc[f"{e}_strip", z] = m.fs.membrane_module.strip_phase.properties[
# #                 0, z
# #             ].conc_mol_comp[e]()
# #             df.loc[f"pH_feed", z] = -log10(
# #                 m.fs.membrane_module.feed_phase.properties[0, z].conc_mol_comp[e]()
# #             )
# #             df.loc[f"pH_strip", z] = -log10(
# #                 m.fs.membrane_module.strip_phase.properties[0, z].conc_mol_comp[e]()
# #             )
# #             df.loc[f"model_pH_feed", z] = m.fs.membrane_module.feed_phase.properties[
# #                 0, z
# #             ].pH_phase()
# #             df.loc[f"model_pH_strip", z] = m.fs.membrane_module.strip_phase.properties[
# #                 0, z
# #             ].pH_phase()
# #         elif e == "H2O":
# #             pass
# #         else:
# #             df.loc[f"{e}_feed", z] = m.fs.membrane_module.feed_phase.properties[
# #                 0, z
# #             ].conc_mol_comp[e]()
# #             df.loc[f"{e}_strip", z] = m.fs.membrane_module.strip_phase.properties[
# #                 0, z
# #             ].conc_mol_comp[e]()

# # df.to_excel("msx_netl_data_steady_profile.xlsx")
