from pyomo.environ import (
    ConcreteModel,
    units,
    Suffix,
    TransformationFactory,
    Var,
    minimize,
    log10,
)
from pyomo.network import Arc, SequentialDecomposition
from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    ControlVolume0DBlock,
)
from idaes.core.util import to_json, from_json, StoreSpec
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.initialization import (
    SingleControlVolumeUnitInitializer,
    BlockTriangularizationInitializer,
)

from prommis.solvent_extraction.membrane_module_property_package import (
    MembraneSXModuleParameters,
)
from prommis.solvent_extraction.non_reactive_tank import NonReactiveTank
from prommis.solvent_extraction.membrane_channel_property_package import (
    MembraneSXChannelParameters,
)
from pyomo.dae.flatten import flatten_dae_components
import matplotlib.pyplot as plt
from prommis.solvent_extraction.membrane_solvent_extraction import (
    MembraneSolventExtraction,
    MembraneSolventExtractionInitializer,
)
import pandas as pd
from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from idaes.core.util import DiagnosticsToolbox

m = ConcreteModel()

time_break_points = [0, 10, 19, 28, 37, 46, 55, 64, 73, 82, 91, 97, 101]
time_fe = 3
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

# for i in range(len(time_break_points)):
#     if i == 0:
#         time_set.append(time_break_points[i])
#     else:
#         mid_point = int((time_break_points[i] + time_break_points[i - 1]) / 2)
#         time_set.append(mid_point)
#         time_set.append(time_break_points[i])

# assert 1 == 2

m.fs = FlowsheetBlock(dynamic=True, time_set=time_set, time_units=units.minute)

m.fs.mem_prop = MembraneSXModuleParameters()
m.fs.mem_channel = MembraneSXChannelParameters()

m.fs.mem_prop.extractant_dosage = 10

m.fs.membrane_module = MembraneSolventExtraction(
    dynamic=False,
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
    finite_elements=25,
    transformation_method="dae.finite_difference",
    transformation_scheme="BACKWARD",
    # collocation_points=2,
    tube_inner_radius=100 * units.micrometer,
    tube_outer_radius=150 * units.micrometer,
    shell_radius=(67 / 2) * units.mm,
    number_of_tubes=6300,
)

m.fs.feed_tank = NonReactiveTank(
    material_balance_type=MaterialBalanceType.componentTotal,
    dynamic=True,
    has_holdup=True,
    property_package=m.fs.mem_channel,
)

m.fs.strip_tank = NonReactiveTank(
    material_balance_type=MaterialBalanceType.componentTotal,
    dynamic=True,
    has_holdup=True,
    property_package=m.fs.mem_channel,
)

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


wts = StoreSpec.value()
from_json(m, fname="membrane_solvent_extraction_5M_HCl.json", wts=wts)
copy_first_steady_state(m)

# Arcs connection

m.fs.feed_to_module = Arc(
    source=m.fs.feed_tank.outlet, destination=m.fs.membrane_module.feed_phase_inlet
)

m.fs.strip_to_module = Arc(
    source=m.fs.strip_tank.outlet, destination=m.fs.membrane_module.strip_phase_inlet
)

m.fs.module_to_strip = Arc(
    source=m.fs.membrane_module.strip_phase_outlet,
    destination=m.fs.strip_tank.inlet,
)

TransformationFactory("network.expand_arcs").apply_to(m)


# Define membrane module parameters

m.fs.membrane_module.module_length.fix(0.254 * 2 * units.m)
m.fs.membrane_module.eff[:, :, "Al"].fix(3.246e-4)
m.fs.membrane_module.eff[:, :, "Ca"].fix(1.765e-3)
m.fs.membrane_module.eff[:, :, "Ce"].fix(0.0305)
m.fs.membrane_module.eff[:, :, "Dy"].fix(1.156e-3)
m.fs.membrane_module.eff[:, :, "Fe"].fix(7.773e-4)
m.fs.membrane_module.eff[:, :, "Gd"].fix(0.013)
m.fs.membrane_module.eff[:, :, "La"].fix(0.032)
m.fs.membrane_module.eff[:, :, "Nd"].fix(0.059)
m.fs.membrane_module.eff[:, :, "Pr"].fix(0.105)
m.fs.membrane_module.eff[:, :, "Sm"].fix(0.15)
m.fs.membrane_module.eff[:, :, "Y"].fix(3e-6)
m.fs.membrane_module.eff[:, :, "Sc"].fix(1.5e-5)

element_to_optimize = ["Pr", "La", "Y", "Ce", "Sm", "Gd", "Nd", "Dy"]

m.a = Var(
    m.fs.mem_prop.component_list, ["c", "l", "q"], initialize=0.1, bounds=(0, None)
)

for e in element_to_optimize:
    m.fs.membrane_module.eff[:, "feed", e].unfix()


@m.Constraint(m.fs.time, element_to_optimize)
def eff_feed(m, t, e):
    v = m.fs.membrane_module.feed_phase.properties[t, 0].flow_vol
    return (
        m.fs.membrane_module.eff[t, "feed", e]
        == m.a[e, "c"] + m.a[e, "l"] * v + m.a[e, "q"] * v**2
    )


# Define feed tank inlet conditions

feed_flow_rate_values = {
    (0, 10): 45,
    (10, 37): 15,
    (37, 64): 30,
    (64, 91): 45,
    (91, 97): 75,
    (97, 101): 150,
}

for t in m.fs.time:
    if t == 0:
        m.fs.feed_tank.inlet.flow_vol[t].fix(
            feed_flow_rate_values[(0, 10)] * units.mL / units.min
        )
    else:
        for k, v in feed_flow_rate_values.items():
            if k[0] < t <= k[1]:
                m.fs.feed_tank.inlet.flow_vol[t].fix(v * units.mL / units.min)

m.fs.feed_tank.inlet.conc_mass_comp[:, "Al"].fix(1351773 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Ca"].fix(6276510 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Fe"].fix(533032 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Sc"].fix(51.4 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "La"].fix(995 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Ce"].fix(1986 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Pr"].fix(275 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Nd"].fix(1104 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Sm"].fix(230 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Gd"].fix(231 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Dy"].fix(199 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "Y"].fix(988 * units.microgram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "H"].fix(10**-2 * 1 * units.gram / units.L)
m.fs.feed_tank.inlet.conc_mass_comp[:, "H2O"].fix(1e6 * units.mg / units.L)
m.fs.feed_tank.control_volume.volume[:].fix(4 * units.L)

# Strip phase inlet conditions

strip_flow_rate_values = {
    (0, 10): 150,
    (10, 19): 45,
    (19, 28): 75,
    (28, 37): 150,
    (37, 46): 45,
    (46, 55): 75,
    (55, 64): 150,
    (64, 73): 45,
    (73, 82): 75,
    (82, 101): 150,
}

for t in m.fs.time:
    if t == 0:
        m.fs.strip_tank.inlet.flow_vol[t].fix(150 * units.mL / units.min)
    else:
        for k, v in strip_flow_rate_values.items():
            if k[0] < t <= k[1]:
                m.fs.strip_tank.inlet.flow_vol[t].fix(v * units.mL / units.min)

m.fs.module_to_strip_expanded.flow_vol_equality.deactivate()

m.fs.strip_tank.control_volume.volume[:].fix(1 * units.L)

# fix conditions at t = 0

# feed tank at t=0
m.fs.feed_tank.control_volume.properties_out[0].flow_vol.fix()
m.fs.strip_tank.control_volume.properties_out[0].flow_vol.fix()


for e in m.fs.mem_channel.component_list:
    if e not in ["H2O"]:
        m.fs.feed_tank.control_volume.properties_out[0].conc_mass_comp[e].fix()
        m.fs.strip_tank.control_volume.properties_out[0].conc_mass_comp[e].fix()


strip_data = pd.read_excel("1M HCl experimental data.xlsx", sheet_name="Sheet3")
strip_data = strip_data.set_index(strip_data.columns[0])

# eff_corr = pd.read_excel("feed_correlation_function.xlsx")
# eff_corr = eff_corr.set_index(eff_corr.columns[0])


@m.Objective(sense=minimize)
def strip_residual(m):
    return (
        (
            sum(
                (
                    m.fs.strip_tank.control_volume.properties_out[t].conc_mass_comp[e]
                    - strip_data[e][t] * 1e-3
                )
                ** 2
                for t in time_break_points
                for e in element_to_optimize
            )
        )
        / (len(time_break_points) * len(element_to_optimize))
    ) ** 0.5


solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 4000
results = solver.solve(m, tee=True)


# to_json(m, fname="membrane_solvent_extraction.json")

for e in element_to_optimize:
    plt.plot(
        m.fs.time, m.fs.strip_tank.control_volume.properties_out[:].conc_mass_comp[e]()
    )
    plt.scatter(time_break_points, [strip_data[e][t] * 1e-3 for t in time_break_points])
    plt.xlabel("Time (min)")
    plt.ylabel(f"{e} strip conc (mg/L)")
    plt.title(f"{e} strip tank conc profile")
    plt.legend(["model", "experiment"])
    plt.show()

eff_corr_5M = pd.DataFrame(columns=["5 M HCl"])
for e in element_to_optimize:
    for s in ["c", "l", "q"]:
        eff_corr_5M.loc[f"{e}_{s}", "5 M HCl"] = m.a[e, s]()

with pd.ExcelWriter(
    "feed_phi_parameters.xlsx", mode="a", engine="openpyxl", if_sheet_exists="replace"
) as writer:
    eff_corr_5M.to_excel(writer, sheet_name="5 M HCl alt", index=True)

flow_alpha = pd.DataFrame(columns=["v", "alpha"], index=[t for t in m.fs.time])

for t in m.fs.time:
    flow_alpha.loc[t, "v"] = m.fs.membrane_module.feed_phase.properties[t, 0].flow_vol()
    flow_alpha.loc[t, "alpha"] = m.fs.membrane_module.eff[t, "feed", "Pr"]()

with pd.ExcelWriter(
    "feed_phi_parameters.xlsx", mode="a", engine="openpyxl", if_sheet_exists="replace"
) as writer:
    flow_alpha.to_excel(writer, sheet_name="5 M HCl alpha values alt", index=True)

fig, ax = plt.subplots(2, 2, figsize=(10, 7))
fig.suptitle("Pr 5 M HCl strip tank profile")
ax[0, 0].step(m.fs.time, m.fs.feed_tank.inlet.flow_vol[:]())
ax[0, 0].set_xlabel("Time, mins")
ax[0, 0].set_ylabel("Feed inlet flowrate, L/hr")
ax[0, 0].set_title("Change in feed inlet flowrate")
ax[0, 1].step(m.fs.time, m.fs.strip_tank.inlet.flow_vol[:]())
ax[0, 1].set_xlabel("Time, mins")
ax[0, 1].set_ylabel("Strip inlet flowrate, L/hr")
ax[0, 1].set_title("Change in strip inlet flowrate")
ax[1, 0].step(m.fs.time, m.fs.feed_tank.inlet.conc_mass_comp[:, "Pr"]())
ax[1, 0].set_xlabel("Time, mins")
ax[1, 0].set_ylabel("Pr feed inlet concentration, mg/L")
ax[1, 0].set_title("Change in Pr feed inlet concentration")
ax[1, 1].plot(
    m.fs.time, m.fs.strip_tank.control_volume.properties_out[:].conc_mass_comp["Pr"]()
)
ax[1, 1].scatter(
    time_break_points, [strip_data["Pr"][t] * 1e-3 for t in time_break_points]
)
ax[1, 1].set_xlabel("Time, mins")
ax[1, 1].set_ylabel("Pr strip tank concentration (mg/L)")
ax[1, 1].set_title("Pr strip tank concentration profile")
ax[1, 1].legend(["model", "experiment"])
plt.tight_layout()
