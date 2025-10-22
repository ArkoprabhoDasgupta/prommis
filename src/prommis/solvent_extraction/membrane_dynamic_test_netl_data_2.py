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
)

from idaes.core.util import DiagnosticsToolbox

m = ConcreteModel()

time_break_points = [0, 12, 22, 32, 42, 52, 63, 72, 82, 92, 102, 117, 122]
time_fe = 2
time_set = []

for i in range(len(time_break_points)):
    if i == 0:
        time_set.append(time_break_points[i])
    else:
        if (time_break_points[i] - time_break_points[i - 1]) >= time_fe:
            number_of_points = int(
                (time_break_points[i] - time_break_points[i - 1]) / time_fe
            )
            if (time_break_points[i] - time_break_points[i - 1]) % time_fe == 0:
                for j in range(1, number_of_points):
                    time_set.append(time_break_points[i - 1] + j * time_fe)
            else:
                for j in range(1, number_of_points + 1):
                    time_set.append(time_break_points[i - 1] + j * time_fe)
            time_set.append(time_break_points[i])
        else:
            time_set.append(time_break_points[i])


m.fs = FlowsheetBlock(dynamic=True, time_set=time_set, time_units=units.minute)

m.fs.mem_prop = MembraneSXModuleParameters()
m.fs.mem_channel = MembraneSXChannelParameters()
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
    finite_elements=2,
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
m.discretizer.apply_to(m, nfe=len(time_set) - 1, wrt=m.fs.time, scheme="BACKWARD")


def copy_first_steady_state(m):
    """
    Function that propagates initial steady state guess to future time points.
    This function is used to initialize all the time discrete variables to the
    initial steady state value.
    """
    regular_vars, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
    # Copy initial conditions forward
    for var in time_vars:
        for t in m.fs.time:
            if t == m.fs.time.first():
                continue
            else:
                var[t].value = var[m.fs.time.first()].value


# from_json(m, fname="membrane_solvent_extraction_dataset_2.json")
# copy_first_steady_state(m)

feed_flow_rate_values = {
    (0, 42): 15,
    (42, 72): 30,
    (72, 102): 45,
    (102, 112): 75,
    (112, 122): 150,
}

feed_concentration_values = {
    (0, 52): {
        "Al": 1338000,
        "Ca": 5720000,
        "Fe": 508344,
        "Y": 1111,
        "La": 1218,
        "Ce": 2506,
        "Pr": 331,
        "Nd": 1332,
        "Sm": 273,
        "Gd": 288,
        "Dy": 239,
        "Sc": 16.6,
    },
    (52, 82): {
        "Al": 1371000,
        "Ca": 5760000,
        "Fe": 468217,
        "Y": 1044,
        "La": 1185,
        "Ce": 2434,
        "Pr": 320,
        "Nd": 1309,
        "Sm": 269,
        "Gd": 274,
        "Dy": 226,
        "Sc": 16.8,
    },
    (82, 102): {
        "Al": 1275000,
        "Ca": 5562000,
        "Fe": 427389,
        "Y": 1017,
        "La": 1136,
        "Ce": 2323,
        "Pr": 306,
        "Nd": 1248,
        "Sm": 258,
        "Gd": 272,
        "Dy": 223,
        "Sc": 21.3,
    },
    (102, 122): {
        "Al": 1334000,
        "Ca": 5897000,
        "Fe": 448437,
        "Y": 1037,
        "La": 1145,
        "Ce": 2393,
        "Pr": 313,
        "Nd": 1268,
        "Sm": 261,
        "Gd": 274,
        "Dy": 221,
        "Sc": 18.9,
    },
}

# Feed phase inlet conditions

for t in m.fs.time:
    if t == 0:
        m.fs.membrane_module.feed_phase_inlet.flow_vol[t].fix(
            feed_flow_rate_values[(0, 42)] * units.mL / units.min
        )
    else:
        for k, v in feed_flow_rate_values.items():
            if k[0] < t <= k[1]:
                m.fs.membrane_module.feed_phase_inlet.flow_vol[t].fix(
                    v * units.mL / units.min
                )

    for k in feed_concentration_values.keys():
        for e, v in feed_concentration_values[k].items():
            if t == 0:
                m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[t, e].fix(
                    feed_concentration_values[(0, 52)][e] * units.microgram / units.L
                )
            elif k[0] < t <= k[1]:
                m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[t, e].fix(
                    v * units.microgram / units.L
                )

m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "H"].fix(
    10**-2.06 * 1 * units.gram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "H2O"].fix(
    1e6 * units.mg / units.L
)

# Strip phase inlet conditions

strip_flow_rate_values = {
    (0, 12): 45,
    (12, 22): 75,
    (22, 42): 150,
    (42, 52): 45,
    (52, 63): 75,
    (63, 72): 150,
    (72, 82): 45,
    (82, 92): 75,
    (92, 122): 150,
}

for t in m.fs.time:
    if t == 0:
        m.fs.membrane_module.strip_phase_inlet.flow_vol[t].fix(
            45 * units.mL / units.min
        )
    else:
        for k, v in strip_flow_rate_values.items():
            if k[0] < t <= k[1]:
                m.fs.membrane_module.strip_phase_inlet.flow_vol[t].fix(
                    v * units.mL / units.min
                )


m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Ca"].fix(
    578.1 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Al"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Fe"].fix(
    92 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "La"].fix(
    0.211 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Ce"].fix(
    0.61 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Pr"].fix(
    0.029 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Nd"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Sm"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Gd"].fix(
    0.013 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Dy"].fix(
    0.046 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Y"].fix(
    0.102 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "Sc"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "H"].fix(
    1 * units.gram / units.L
)

m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[:, "H2O"].fix(
    1e6 * units.mg / units.L
)

m.fs.membrane_module.eff["Al"].fix(3.246e-4)
m.fs.membrane_module.eff["Ca"].fix(1.765e-3)
m.fs.membrane_module.eff["Ce"].fix(0.0305)
m.fs.membrane_module.eff["Dy"].fix(1.156e-3)
m.fs.membrane_module.eff["Fe"].fix(7.773e-4)
m.fs.membrane_module.eff["Gd"].fix(0.013)
m.fs.membrane_module.eff["La"].fix(0.032)
m.fs.membrane_module.eff["Nd"].fix(0.059)
m.fs.membrane_module.eff["Pr"].fix(0.105)
m.fs.membrane_module.eff["Sm"].fix(0.15)
m.fs.membrane_module.eff["Y"].fix(3e-6)
m.fs.membrane_module.eff["Sc"].fix(1.5e-5)

# for e in m.fs.mem_channel.component_list:
#     if e not in ["H2O"]:
#         m.fs.membrane_module.feed_phase.properties[0.0, :].conc_mass_comp[e].fix()
#         m.fs.membrane_module.feed_phase.properties[0.0, :].flow_vol.fix()

#         m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp[e].fix()
#         m.fs.membrane_module.strip_phase.properties[0.0, :].flow_vol.fix()


m.fs.membrane_module.strip_phase.properties[0.0, :].flow_vol.fix(
    45 * units.mL / units.min
)

m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Ca"].fix(
    578.1 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Al"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Fe"].fix(
    92 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["La"].fix(
    0.211 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Ce"].fix(
    0.61 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Pr"].fix(
    0.029 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Nd"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Sm"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Gd"].fix(
    0.013 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Dy"].fix(
    0.046 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Y"].fix(
    0.102 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["Sc"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["H"].fix(
    1 * units.gram / units.L
)

# m.fs.membrane_module.strip_phase.properties[0.0, :].conc_mass_comp["H2O"].fix(
#     1e6 * units.mg / units.L
# )


for k in feed_concentration_values.keys():
    for e, v in feed_concentration_values[k].items():
        m.fs.membrane_module.feed_phase.properties[0.0, :].conc_mass_comp[e].fix(
            feed_concentration_values[(0, 52)][e] * units.microgram / units.L
        )

# m.fs.membrane_module.feed_phase.properties[0.0, :].conc_mass_comp["H2O"].fix(
#     1e6 * units.mg / units.L
# )
m.fs.membrane_module.feed_phase.properties[0, :].conc_mass_comp["H"].fix(
    10**-2.06 * 1 * units.gram / units.L
)
m.fs.membrane_module.feed_phase.properties[0.0, :].flow_vol.fix(
    15 * units.mL / units.min
)

print("Degrees of freedom:", dof(m))

# m.scaling_factor = Suffix(direction=Suffix.EXPORT)

# for t in m.fs.time:
#     for z in m.fs.membrane_module.feed_phase.length_domain:
#         set_scaling_factor(
#             m.fs.membrane_module.feed_phase.properties[t, z].conc_mass_comp["H"],
#             1e2,
#         )
#         set_scaling_factor(
#             m.fs.membrane_module.feed_phase.properties[t, z].pH_constraint,
#             1e2,
#         )
#         set_scaling_factor(
#             m.fs.membrane_module.feed_phase.material_flow_linking_constraints[
#                 t, z, "liquid", "H"
#             ],
#             1e2,
#         )
#         for e in m.fs.mem_prop.component_list:
#             set_scaling_factor(
#                 m.fs.membrane_module.feed_phase.properties[t, z].conc_mol_comp[e],
#                 1e2,
#             )
#             set_scaling_factor(
#                 m.fs.membrane_module.feed_phase.properties[t, z].conc_mass_comp[e],
#                 1,
#             )
#             set_scaling_factor(
#                 m.fs.membrane_module.strip_phase.properties[t, z].conc_mol_comp[e],
#                 1e2,
#             )
#             set_scaling_factor(
#                 m.fs.membrane_module.strip_phase.properties[t, z].conc_mass_comp[e],
#                 1,
#             )

#             for r in m.fs.membrane_module.r:
#                 set_scaling_factor(
#                     m.fs.membrane_module.conc_mol_membrane_comp[t, z, r, e],
#                     1e2,
#                 )


# scaling = TransformationFactory("core.scale_model")
# scaled_model = scaling.create_using(m, rename=False)

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 10000
results = solver.solve(m, tee=True)
# results = solver.solve(scaled_model, tee=True)

# scaling.propagate_solution(scaled_model, m)

strip_percentage_recovery = {}

for e in m.fs.mem_prop.component_list:
    strip_percentage_recovery[e] = [
        (
            (
                (
                    m.fs.membrane_module.strip_phase.properties[0, z].conc_mass_comp[
                        e
                    ]()
                    * m.fs.membrane_module.strip_phase.properties[0, z].flow_vol()
                )
                - (
                    m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, e]()
                    * m.fs.membrane_module.strip_phase_inlet.flow_vol[0]()
                )
            )
            / (
                m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, e]()
                * m.fs.membrane_module.feed_phase_inlet.flow_vol[0]()
            )
        )
        * 100
        for z in m.fs.membrane_module.strip_phase.length_domain
    ]

strip_outlet_recovery = {}

for e in m.fs.mem_prop.component_list:
    strip_outlet_recovery[e] = [
        (
            (
                (
                    m.fs.membrane_module.strip_phase_outlet.conc_mass_comp[t, e]()
                    * m.fs.membrane_module.strip_phase_outlet.flow_vol[t]()
                )
                - (
                    m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[t, e]()
                    * m.fs.membrane_module.strip_phase_inlet.flow_vol[t]()
                )
            )
            / (
                m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[t, e]()
                * m.fs.membrane_module.feed_phase_inlet.flow_vol[t]()
            )
        )
        * 100
        for t in m.fs.time
    ]

import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 2, figsize=(10, 7))
ax[0, 0].plot(m.fs.time, m.fs.membrane_module.feed_phase_inlet.flow_vol[:]())
ax[0, 0].set_xlabel("Time, mins")
ax[0, 0].set_ylabel("Feed inlet flowrate, L/hr")
ax[0, 0].set_title("Change in feed inlet flowrate")
ax[0, 1].plot(m.fs.time, m.fs.membrane_module.strip_phase_inlet.flow_vol[:]())
ax[0, 1].set_xlabel("Time, mins")
ax[0, 1].set_ylabel("Strip inlet flowrate, L/hr")
ax[0, 1].set_title("Change in strip inlet flowrate")
ax[1, 0].plot(
    m.fs.time, m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[:, "Pr"]()
)
ax[1, 0].set_xlabel("Time, mins")
ax[1, 0].set_ylabel("Pr feed inlet concentration, mg/L")
ax[1, 0].set_title("Change in Pr feed inlet concentration")
ax[1, 1].plot(m.fs.time, strip_outlet_recovery["Pr"])
ax[1, 1].set_xlabel("Time, mins")
ax[1, 1].set_ylabel("Pr strip outlet recovery %")
ax[1, 1].set_title("Change in Pr strip outlet recovery %")
plt.tight_layout()
