from pyomo.environ import (
    ConcreteModel,
    units,
    Suffix,
    TransformationFactory,
    Var,
    minimize,
)

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.util import to_json
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.util.scaling import set_scaling_factor

from prommis.solvent_extraction.membrane_module_property_package import (
    MembraneSXModuleParameters,
)

from prommis.solvent_extraction.membrane_channel_property_package import (
    MembraneSXChannelParameters,
)

# from prommis.solvent_extraction.membrane_channel_property_package_mass_flow import (
#     MembraneSXChannelParameters,
# )
from prommis.solvent_extraction.membrane_solvent_extraction import (
    MembraneSolventExtraction,
    MembraneSolventExtractionInitializer,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from idaes.core.util import DiagnosticsToolbox

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

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
    finite_elements=25,
    transformation_method="dae.finite_difference",
    transformation_scheme="BACKWARD",
    # collocation_points=2,
    tube_inner_radius=100 * units.micrometer,
    tube_outer_radius=150 * units.micrometer,
    shell_radius=(67 / 2) * units.mm,
    number_of_tubes=6300,
)

m.fs.membrane_module.module_length.fix(0.254 * 2 * units.m)

# Feed phase inlet conditions

m.fs.membrane_module.feed_phase_inlet.flow_vol.fix(17 * units.mL / units.min)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Ca"].fix(
    1535976 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Al"].fix(
    310800 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Fe"].fix(
    131300 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "La"].fix(
    273 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Ce"].fix(
    573 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Pr"].fix(
    70 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Nd"].fix(
    315 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Sm"].fix(
    65.9 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Gd"].fix(
    67.5 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Dy"].fix(
    52.4 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Y"].fix(
    231 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Sc"].fix(
    20 * units.microgram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H"].fix(
    10**-2.06 * 1 * units.gram / units.L
)
# m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H"].fix(
#     0.01 * units.gram / units.L
# )
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H2O"].fix(1e6)

# Strip phase inlet conditions

m.fs.membrane_module.strip_phase_inlet.flow_vol.fix(54 * units.mL / units.min)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Ca"].fix(
    4081 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Al"].fix(
    1e-7 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Fe"].fix(
    515 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "La"].fix(
    1.89 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Ce"].fix(
    3.19 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Pr"].fix(
    0.398 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Nd"].fix(
    1.24 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Sm"].fix(
    0.15 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Gd"].fix(
    0.157 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Dy"].fix(
    0.119 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Y"].fix(
    0.6 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Sc"].fix(
    0.649 * units.microgram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "H"].fix(
    3 * 1 * units.gram / units.L
)
# m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "H"].fix(
#     0.5 * units.gram / units.L
# )
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "H2O"].fix(1e6)

# m.fs.membrane_module.eff[:].fix(1)

m.fs.membrane_module.eff["Al"].fix(3.246e-4)
m.fs.membrane_module.eff["Ca"].fix(1.765e-3)
m.fs.membrane_module.eff["Ce"].fix(0.02709)
m.fs.membrane_module.eff["Dy"].fix(1.156e-4)
m.fs.membrane_module.eff["Fe"].fix(7.773e-4)
m.fs.membrane_module.eff["Gd"].fix(0.00998)
m.fs.membrane_module.eff["La"].fix(0.0291)
m.fs.membrane_module.eff["Nd"].fix(0.05226)
m.fs.membrane_module.eff["Pr"].fix(0.0993)
m.fs.membrane_module.eff["Sm"].fix(0.124)
m.fs.membrane_module.eff["Y"].fix(7.42e-7)

m.fs.membrane_module.eff["Sc"].fix(0.4)

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

# initializer = MembraneSolventExtractionInitializer()
# initializer.initialize(scaled_model.fs.membrane_module)
# # initializer.initialize(m.fs.membrane_module)

print("Degrees of freedom:", dof(m))

# m.final_extraction_value = Var(m.fs.mem_prop.component_list, initialize=100)


# @m.Constraint(m.fs.mem_prop.component_list)
# def final_ext_func(m, e):
#     return (
#         m.final_extraction_value[e]
#         == (
#             1
#             - (
#                 m.fs.membrane_module.feed_phase.properties[0, 1].conc_mass_comp[e]
#                 * m.fs.membrane_module.feed_phase.properties[0, 1].flow_vol
#             )
#             / (
#                 m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, e]
#                 * m.fs.membrane_module.feed_phase_inlet.flow_vol[0]
#             )
#         )
#         * 100
#     )


final_ext_value = {
    "Al": 0.3,
    "Ca": 2,
    "Fe": 2,
    "Sc": 71,
    "La": 53.3,
    "Ce": 82.6,
    "Pr": 87.8,
    "Nd": 89.5,
    "Sm": 97.5,
    "Gd": 96.1,
    "Dy": 99,
    "Y": 99.3,
}


print("Degrees of freedom:", dof(m))

# @m.Objective(sense=minimize)
# def objective_function(m):

#     return (
#         sum(
#             (m.final_extraction_value[e] - final_ext_value[e]) ** 2
#             for e in m.fs.mem_prop.component_list
#         )
#     ) ** 0.5


# @m.Objective(sense=minimize)
# def objective_function(m):

#     return (((m.final_extraction_value["Sc"] - final_ext_value["Sc"]) ** 2)) ** 0.5

# initializer = MembraneSolventExtractionInitializer()
# initializer.initialize(m.fs.membrane_module)

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 2000
results = solver.solve(m, tee=True)
# results = solver.solve(scaled_model, tee=True)

# scaling.propagate_solution(scaled_model, m)

feed_percentage_recovery = {}
strip_percentage_recovery = {}

for e in m.fs.mem_prop.component_list:
    feed_percentage_recovery[e] = [
        (
            1
            - (
                m.fs.membrane_module.feed_phase.properties[0, z].conc_mass_comp[e]()
                * m.fs.membrane_module.feed_phase.properties[0, z].flow_vol()
            )
            / (
                m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, e]()
                * m.fs.membrane_module.feed_phase_inlet.flow_vol[0]()
            )
        )
        * 100
        for z in m.fs.membrane_module.feed_phase.length_domain
    ]

for e in m.fs.mem_prop.component_list:
    strip_percentage_recovery[e] = [
        (
            (
                m.fs.membrane_module.strip_phase.properties[0, z].conc_mass_comp[e]()
                * m.fs.membrane_module.strip_phase.properties[0, z].flow_vol()
            )
            / (
                m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, e]()
                * m.fs.membrane_module.strip_phase_inlet.flow_vol[0]()
            )
            - 1
        )
        * 100
        for z in m.fs.membrane_module.strip_phase.length_domain
    ]
