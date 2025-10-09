import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    TransformationFactory,
    units,
    value,
    Var,
)
from pyomo.common.collections import ComponentSet
from pyomo.util.subsystems import (
    TemporarySubsystemManager,
    create_subsystem_block,
)

from pyomo.contrib.incidence_analysis import IncidenceGraphInterface

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.core.solvers import get_solver, petsc
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.model_serializer as ms

from idaes.models.unit_models.cstr import CSTR
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)

from prommis.solvent_extraction.dae_utils import (
    generate_diff_deriv_disc_components_along_set,
    DAEInterface,
)


def pause():
    programPause = input("Press the <ENTER> key to continue...")


def fix_initial_conditions(model, time_set, t0=None):
    if t0 is None:
        t0 = time_set.first()
    for state, _, _ in generate_diff_deriv_disc_components_along_set(model, time_set):
        state[t0].fix()


def create_model(time_set=None):
    # CSTR model setup copy and pasted from test_cstr.py (the original authors are Andrew and Vibhav)

    m = ConcreteModel()
    if time_set is None:
        m.fs = FlowsheetBlock(dynamic=False)
    else:
        m.fs = FlowsheetBlock(dynamic=True, time_set=time_set, time_units=units.s)

    m.fs.properties = SaponificationParameterBlock()
    m.fs.reactions = SaponificationReactionParameterBlock(
        property_package=m.fs.properties
    )

    m.fs.unit = CSTR(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_equilibrium_reactions=False,
        has_heat_transfer=True,
        has_heat_of_reaction=True,
        has_pressure_change=True,
    )

    m.fs.unit.inlet.flow_vol.fix(1.0e-03)
    m.fs.unit.inlet.conc_mol_comp[:, "H2O"].fix(55388.0)
    m.fs.unit.inlet.conc_mol_comp[:, "NaOH"].fix(100.0)
    m.fs.unit.inlet.conc_mol_comp[:, "EthylAcetate"].fix(100.0)
    m.fs.unit.inlet.conc_mol_comp[:, "SodiumAcetate"].fix(1e-8)
    m.fs.unit.inlet.conc_mol_comp[:, "Ethanol"].fix(1e-8)

    m.fs.unit.inlet.temperature.fix(303.15)
    m.fs.unit.inlet.pressure.fix(101325.0)

    m.fs.unit.volume.fix(1.5e-03)
    m.fs.unit.heat_duty.fix(0)
    m.fs.unit.deltaP.fix(0)

    return m


if __name__ == "__main__":
    index_reduction = True

    m_steady = create_model()

    assert degrees_of_freedom(m_steady) == 0
    init = m_steady.fs.unit.default_initializer()
    init.initialize(m_steady.fs.unit)
    solver_obj = get_solver("ipopt_v2")

    results = solver_obj.solve(m_steady, tee=True)
    assert_optimal_termination(results)

    initial_condition = ms.to_json(m_steady, return_dict=True)

    if not index_reduction:
        m_dynamic = create_model(time_set=[60 * i for i in range(11)])

        ms.from_json(m_dynamic, initial_condition)

        TransformationFactory("dae.finite_difference").apply_to(
            m_dynamic,
            nfe=len(m_dynamic.fs.time) - 1,
            wrt=m_dynamic.fs.time,
            scheme="BACKWARD",
        )

        print(degrees_of_freedom(m_dynamic))
        fix_initial_conditions(m_dynamic, m_dynamic.fs.time)
        print(degrees_of_freedom(m_dynamic))

        dae_dynamic = DAEInterface(m_dynamic, m_dynamic.fs.time)

        # alg_vars, alg_cons = dae.get_algebraic_subsystem_at_time(t=60)
        alg_vars, alg_cons = dae_dynamic.get_naive_algebraic_subsystem_at_time(t=60)

        igraph = IncidenceGraphInterface()

        var_dmp, con_dmp = igraph.dulmage_mendelsohn(
            variables=alg_vars, constraints=alg_cons
        )

        print("Unpaired Variables")
        for var in var_dmp[0]:
            print(var.name)
        print("\n")
        print("Underconstrained Variables")
        for var in var_dmp[1]:
            print(var.name)
        print("\n")
        print("Overconstrained Variables")
        for var in var_dmp[2]:
            print(var.name)
        print("\n\n")

        print("Unpaired Constraints")
        for con in con_dmp[0]:
            print(con.name)
        print("\n")
        print("Underconstrained Constraints")
        for con in con_dmp[2]:
            print(con.name)
        print("\n")
        print("Overconstrained Constraints")
        for con in con_dmp[1]:
            print(con.name)
        print("\n\n")

        pause()

        m_dynamic.fs.unit.inlet.flow_vol.fix(1.5e-03)

        petsc.petsc_dae_by_time_element(
            m_dynamic,
            time=m_dynamic.fs.time,
            between=[m_dynamic.fs.time.first(), m_dynamic.fs.time.last()],
            ts_options={
                "--ts_type": "beuler",
                "--ts_dt": 0.1,
                "--ts_monitor": "",  # set initial step to 0.1
                "--ts_save_trajectory": 1,
            },
            symbolic_solver_labels=True,
        )

    else:
        m_reduced = create_model(time_set=[60 * i for i in range(11)])

        ms.from_json(m_reduced, initial_condition)

        TransformationFactory("dae.finite_difference").apply_to(
            m_reduced,
            nfe=len(m_reduced.fs.time) - 1,
            wrt=m_reduced.fs.time,
            scheme="BACKWARD",
        )

        print(degrees_of_freedom(m_reduced))

        m_reduced.fs.unit.control_volume.material_accumulation_disc_eq[
            :, "Liq", "H2O"
        ].deactivate()

        @m_reduced.fs.unit.control_volume.Constraint(m_reduced.fs.time)
        def volume_balance_eqn(b, t):
            return b.material_accumulation[t, "Liq", "H2O"] == 0

        print(degrees_of_freedom(m_reduced))
        fix_initial_conditions(m_reduced, m_reduced.fs.time)
        print(degrees_of_freedom(m_reduced))

        dae_reduced = DAEInterface(m_reduced, m_reduced.fs.time)

        # alg_vars, alg_cons = dae.get_algebraic_subsystem_at_time(t=60)
        alg_vars, alg_cons = dae_reduced.get_naive_algebraic_subsystem_at_time(t=60)

        igraph = IncidenceGraphInterface()

        var_dmp, con_dmp = igraph.dulmage_mendelsohn(
            variables=alg_vars, constraints=alg_cons
        )

        print("Unpaired Variables")
        for var in var_dmp[0]:
            print(var.name)
        print("\n")
        print("Underconstrained Variables")
        for var in var_dmp[1]:
            print(var.name)
        print("\n")
        print("Overconstrained Variables")
        for var in var_dmp[2]:
            print(var.name)
        print("\n\n")

        print("Unpaired Constraints")
        for con in con_dmp[0]:
            print(con.name)
        print("\n")
        print("Underconstrained Constraints")
        for con in con_dmp[2]:
            print(con.name)
        print("\n")
        print("Overconstrained Constraints")
        for con in con_dmp[1]:
            print(con.name)
        print("\n\n")

        pause()

        m_reduced.fs.unit.inlet.flow_vol.fix(1.5e-03)

        petsc.petsc_dae_by_time_element(
            m_reduced,
            time=m_reduced.fs.time,
            between=[m_reduced.fs.time.first(), m_reduced.fs.time.last()],
            ts_options={
                "--ts_type": "beuler",
                "--ts_dt": 0.1,
                "--ts_monitor": "",  # set initial step to 0.1
                "--ts_save_trajectory": 1,
            },
            symbolic_solver_labels=True,
        )
