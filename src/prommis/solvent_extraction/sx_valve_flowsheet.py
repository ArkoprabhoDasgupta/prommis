from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    units,
    TransformationFactory,
    Var,
    value,
    log10,
    Suffix,
    Constraint,
)
from pyomo.dae.flatten import flatten_dae_components
from pyomo.network import Arc

import numpy as np

import matplotlib.pyplot as plt

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
import idaes.core.solvers.petsc as petsc
from idaes.core.util.initialization import initialize_by_time_element
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util import from_json, DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.scaling import (
    set_scaling_factor,
    CustomScalerBase,
    report_scaling_factors,
)
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Valve
from idaes.models.control.controller import (
    PIDController,
    ControllerType,
    ControllerMVBoundType,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.hybrid_solvent_extraction import SolventExtraction
from prommis.solvent_extraction.D_variant_model import D_calculation

m = ConcreteModel()

time_duration = 24

Elements = ["Y", "Ce", "Nd", "Sm", "Gd", "Dy"]

m.fs = FlowsheetBlock(dynamic=True, time_set=[0, time_duration], time_units=units.hour)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()

number_of_stages = 1
stage_number = np.arange(1, number_of_stages + 1)

m.fs.solex_1 = SolventExtraction(
    number_of_finite_elements=number_of_stages,
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
)

m.fs.solex_2 = SolventExtraction(
    number_of_finite_elements=number_of_stages,
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
)


def _valve_pressure_flow_cb(b):

    b.Cv = Var(initialize=1e-2, units=units.dm**2)
    b.Cv.fix()

    @b.Constraint(b.flowsheet().time)
    def pressure_flow_equation(b, t):
        # rho_aqueous = units.convert(
        #     b.control_volume.properties_in[t].dens_mass,
        #     to_units=units.kg / (units.m**3),
        # )
        rho_aqueous = units.convert(
            sum(
                b.control_volume.properties_out[t].conc_mass_comp[p]
                for p in m.fs.leach_soln.component_list
            ),
            to_units=units.kg / (units.m**3),
        )
        Po = units.convert(
            b.control_volume.properties_out[t].pressure, to_units=units.Pa
        )
        Pi = units.convert(
            b.control_volume.properties_in[t].pressure, to_units=units.Pa
        )
        vel_head = units.convert(
            (Pi - Po) / rho_aqueous, to_units=(units.dm**2 / units.hr**2)
        )

        F = units.convert(
            b.control_volume.properties_out[t].flow_vol,
            to_units=units.L / units.hr,
        )
        Cv = b.Cv
        fun = b.valve_function[t]
        return F**2 == ((Cv * fun) ** 2) * vel_head
        # return (Cv * F) ** 2 == vel_head * fun**2


m.fs.valve_1 = Valve(
    dynamic=False,
    has_holdup=False,
    material_balance_type=MaterialBalanceType.componentTotal,
    property_package=m.fs.leach_soln,
    pressure_flow_callback=_valve_pressure_flow_cb,
)

m.fs.valve_2 = Valve(
    dynamic=False,
    has_holdup=False,
    material_balance_type=MaterialBalanceType.componentTotal,
    property_package=m.fs.leach_soln,
    pressure_flow_callback=_valve_pressure_flow_cb,
)

m.fs.sx_1_to_v_1 = Arc(
    source=m.fs.solex_1.aqueous_outlet, destination=m.fs.valve_1.inlet
)
m.fs.v_1_to_sx_2 = Arc(
    source=m.fs.valve_1.outlet, destination=m.fs.solex_2.aqueous_inlet
)
m.fs.sx_2_to_v_2 = Arc(
    source=m.fs.solex_2.aqueous_outlet, destination=m.fs.valve_2.inlet
)


@m.Constraint(m.fs.time, stage_number)
def pressure_calculation_1(m, t, s):
    g = 9.8 * (units.m) / units.sec**2
    P_atm = 101325 * units.Pa
    A = 1 * units.m**2

    rho_aq = sum(
        m.fs.solex.mscontactor.aqueous[t, s].conc_mass_comp[p]
        for p in m.fs.leach_soln.component_list
    )
    rho_og = sum(
        m.fs.solex.mscontactor.organic[t, s].conc_mass_comp[p]
        for p in m.fs.prop_o.component_list
    )
    P_aq = units.convert(
        (
            rho_aq
            * g
            * m.fs.solex.mscontactor.volume[s]
            * m.fs.solex.mscontactor.volume_frac_stream[t, s, "aqueous"]
            / A
        ),
        to_units=units.Pa,
    )
    P_org = units.convert(
        (
            rho_og
            * g
            * m.fs.solex.mscontactor.volume[s]
            * m.fs.solex.mscontactor.volume_frac_stream[t, s, "organic"]
            / A
        ),
        to_units=units.Pa,
    )
    hydrostat_head = 10 * units.m
    hydrostat_pressure = units.convert(
        (rho_og * g * hydrostat_head),
        to_units=units.Pa,
    )
    return (
        m.fs.solex.mscontactor.aqueous[t, s].pressure
        == P_aq + P_org + P_atm + hydrostat_pressure
    )
