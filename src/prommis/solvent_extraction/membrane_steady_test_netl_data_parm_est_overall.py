from pyomo.environ import ConcreteModel, units, minimize, value
import numpy as np
from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom as dof

from prommis.solvent_extraction.membrane_module_property_package_parm_est import (
    MembraneSXModuleParameters,
)
from prommis.solvent_extraction.membrane_solvent_extraction import (
    MembraneSolventExtraction,
    MembraneSolventExtractionInitializer,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.mem_prop = MembraneSXModuleParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.mem_prop.extractant_dosage = 10

m.fs.membrane_module = MembraneSolventExtraction(
    feed_phase={
        "property_package": m.fs.leach_soln,
        "flow_direction": FlowDirection.forward,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "energy_balance_type": EnergyBalanceType.none,
    },
    strip_phase={
        "property_package": m.fs.leach_soln,
        "flow_direction": FlowDirection.forward,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "energy_balance_type": EnergyBalanceType.none,
    },
    membrane_phase={
        "property_package": m.fs.mem_prop,
    },
    finite_elements=6,
    transformation_method="dae.finite_difference",
    transformation_scheme="BACKWARD",
    tube_inner_radius=100 * units.micrometer,
    tube_outer_radius=150 * units.micrometer,
    shell_radius=(67 / 2) * units.mm,
    number_of_tubes=6300,
)

m.fs.membrane_module.module_length.fix(0.254 * 2 * units.m)

# Feed phase inlet conditions

m.fs.membrane_module.feed_phase_inlet.flow_vol.fix(15 * units.mL / units.min)
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
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "Cl"].fix(
    0.01 * 35.5 * units.gram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H"].fix(
    10**-2.06 * 1 * units.gram / units.L
)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "HSO4"].fix(1e-5)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "SO4"].fix(1e-5)
m.fs.membrane_module.feed_phase_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.membrane_module.feed_phase_inlet.temperature.fix(303)
m.fs.membrane_module.feed_phase_inlet.pressure.fix(101325)

# Strip phase inlet conditions

m.fs.membrane_module.strip_phase_inlet.flow_vol.fix(10 * units.mL / units.min)
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
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "Cl"].fix(
    3 * 35.5 * units.gram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "H"].fix(
    3 * 1 * units.gram / units.L
)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "HSO4"].fix(1e-5)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "SO4"].fix(1e-5)
m.fs.membrane_module.strip_phase_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.membrane_module.strip_phase_inlet.temperature.fix(303)
m.fs.membrane_module.strip_phase_inlet.pressure.fix(101325)

m.fs.membrane_module.feed_phase.properties[:, :].temperature.fix(303)
m.fs.membrane_module.feed_phase.properties[:, :].pressure.fix(101325)
m.fs.membrane_module.strip_phase.properties[:, :].temperature.fix(303)
m.fs.membrane_module.strip_phase.properties[:, :].pressure.fix(101325)

m.fs.mem_prop.D_coeff["Sc"].fix()


feed_outlet_value = {
    "Al": 263600 * units.microgram / units.L,
    "Ca": 1211976 * units.microgram / units.L,
    "Fe": 99370 * units.microgram / units.L,
    # "Sc": 3.79 * units.microgram / units.L,
    "La": 132 * units.microgram / units.L,
    "Ce": 15.2 * units.microgram / units.L,
    "Pr": 1.05 * units.microgram / units.L,
    "Nd": 3.21 * units.microgram / units.L,
    "Sm": 0.171 * units.microgram / units.L,
    "Gd": 0.351 * units.microgram / units.L,
    "Dy": 0.162 * units.microgram / units.L,
    "Y": 0.727 * units.microgram / units.L,
}

strip_outlet_value = {
    "Al": 838 * units.microgram / units.L,
    "Ca": 27946 * units.microgram / units.L,
    "Fe": 1549 * units.microgram / units.L,
    # "Sc": 3.79 * units.microgram / units.L,
    "La": 88.2 * units.microgram / units.L,
    "Ce": 276 * units.microgram / units.L,
    "Pr": 37.8 * units.microgram / units.L,
    "Nd": 158 * units.microgram / units.L,
    "Sm": 29.5 * units.microgram / units.L,
    "Gd": 29.2 * units.microgram / units.L,
    "Dy": 14.4 * units.microgram / units.L,
    "Y": 35.1 * units.microgram / units.L,
}

optim_comp_list = m.fs.mem_prop.component_list - ["Sc"]


@m.Expression(optim_comp_list)
def feed_deviation(m, element):
    return abs(
        m.fs.membrane_module.feed_phase.properties[0, 1].conc_mass_comp[element]
        - units.convert(feed_outlet_value[element], to_units=units.mg / units.L)
    )


@m.Expression(optim_comp_list)
def strip_deviation(m, element):
    return abs(
        m.fs.membrane_module.strip_phase.properties[0, 1].conc_mass_comp[element]
        - units.convert(strip_outlet_value[element], to_units=units.mg / units.L)
    )


@m.Objective(sense=minimize)
def objective_function(m):

    return (
        sum(value(m.feed_deviation[e]) for e in optim_comp_list)
        + sum(value(m.strip_deviation[e]) for e in optim_comp_list)
    ) / 22


print(dof(m))
solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 100000
results = solver.solve(m, tee=True)
