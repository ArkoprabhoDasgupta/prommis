from pyomo.environ import ConcreteModel, SolverFactory, Constraint, units, log10
from pyomo.environ import TransformationFactory

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.neutralization_tank import NeutralizationTank
from idaes.core.scaling import CustomScalerBase, report_scaling_factors

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.leach_soln = LeachSolutionParameters()

m.fs.neutral = NeutralizationTank(property_package=m.fs.leach_soln)

m.fs.neutral.base_flowrate[0].fix(2)
m.fs.neutral.base_concentration[0].fix(0.2)

m.fs.neutral.inlet.flow_vol.fix(62.01)

m.fs.neutral.inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.neutral.inlet.conc_mass_comp[0, "H"].fix(1.755)
m.fs.neutral.inlet.conc_mass_comp[0, "SO4"].fix(3999.818)
m.fs.neutral.inlet.conc_mass_comp[0, "HSO4"].fix(693.459)
m.fs.neutral.inlet.conc_mass_comp[0, "Al"].fix(422.375)
m.fs.neutral.inlet.conc_mass_comp[0, "Ca"].fix(109.542)
m.fs.neutral.inlet.conc_mass_comp[0, "Fe"].fix(688.266)
m.fs.neutral.inlet.conc_mass_comp[0, "Sc"].fix(0.032)
m.fs.neutral.inlet.conc_mass_comp[0, "Y"].fix(0.124)
m.fs.neutral.inlet.conc_mass_comp[0, "La"].fix(0.986)
m.fs.neutral.inlet.conc_mass_comp[0, "Ce"].fix(2.277)
m.fs.neutral.inlet.conc_mass_comp[0, "Pr"].fix(0.303)
m.fs.neutral.inlet.conc_mass_comp[0, "Nd"].fix(0.946)
m.fs.neutral.inlet.conc_mass_comp[0, "Sm"].fix(0.097)
m.fs.neutral.inlet.conc_mass_comp[0, "Gd"].fix(0.2584)
m.fs.neutral.inlet.conc_mass_comp[0, "Dy"].fix(0.047)
m.fs.neutral.inlet.conc_mass_comp[0, "Cl"].fix(1e-8)


class NeutralizeScale(CustomScalerBase):

    DEFAULT_SCALING_FACTORS = {"flow_vol": 1e-2, "conc_mass_comp": 1e-2}

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        self.scale_variable_by_default(model.flow_vol, overwrite=overwrite)
        for k, v in model.conc_mass_comp.items():
            if k == "H2O":
                self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
            else:
                self.scale_variable_by_default(v, overwrite=False)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        if model.is_property_constructed("h2o_concentration"):
            self.set_constraint_scaling_factor(
                model.h2o_concentration, 1e-4, overwrite=overwrite
            )
        if model.is_property_constructed("molar_concentration_constraint"):
            for v in model.molar_concentration_constraint.values():
                self.set_constraint_scaling_factor(v, 1e-3, overwrite=overwrite)


class NeutralTankScale(CustomScalerBase):

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):

        self.call_submodel_scaler_method(
            model=model,
            submodel="control_volume.properties_in",
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        self.propagate_state_scaling(
            target_state=model.control_volume.properties_out,
            source_state=model.control_volume.properties_in,
            overwrite=overwrite,
        )

        self.call_submodel_scaler_method(
            model=model,
            submodel="control_volume.properties_out",
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        for v in model.base_flowrate.values():
            self.set_variable_scaling_factor(v, 1e-2)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Call scaling methods for sub-models
        self.call_submodel_scaler_method(
            model=model,
            submodel="control_volume.properties_in",
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            model=model,
            submodel="control_volume.properties_out",
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Scale control volume constraints
        for c in model.control_volume.component_data_objects(
            Constraint, descend_into=False
        ):
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )

        # Scale main model constraints
        if hasattr(model, "mass_transfer_constraint"):
            for c in model.mass_transfer_constraint.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )


scaler = NeutralTankScale()
scaler.scale_model(
    m.fs.neutral,
    submodel_scalers={
        "control_volume.properties_in": NeutralizeScale,
        "control_volume.properties_out": NeutralizeScale,
    },
)

report_scaling_factors(m, descend_into=True)

solver = get_solver(
    "ipopt_v2", writer_config={"linear_presolve": True, "scale_model": True}
)
solver.solve(m, tee=True)

print(-log10(m.fs.neutral.control_volume.properties_in[0].conc_mol_comp["H"]()))
print(-log10(m.fs.neutral.control_volume.properties_out[0].conc_mol_comp["H"]()))

sm = TransformationFactory("core.scale_model").create_using(m, rename=False)
dt = DiagnosticsToolbox(sm)
