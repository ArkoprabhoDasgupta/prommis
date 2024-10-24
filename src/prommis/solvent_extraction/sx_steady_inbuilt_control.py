#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

""" 
Demonstration flowsheet for dynamic state solvent extraction loading process
using parameters and data derived from West Kentucky No. 13 coal refuse.

Authors: Arkoprabho Dasgupta

"""

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    units,
    TransformationFactory,
    Var,
    value,
    log10,
    Suffix,
    Constraint
)

import numpy as np

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.util import DiagnosticsToolbox
from idaes.core.scaling import CustomScalerBase, report_scaling_factors, set_scaling_factor
from idaes.core.solvers import get_solver

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.hybrid_solvent_extraction import SolventExtraction
from prommis.solvent_extraction.D_variant_model import D_calculation

"""
Method of building a dynamic solvent extraction model with a specified number of
stages and with two separate property packages for the two inlet streams.
This is a loading operation, so no additional argument has to be specified.

"""

m = ConcreteModel()

m.fs = FlowsheetBlock()

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()

number_of_stages = 3

m.fs.solex = SolventExtraction(
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


stage_number = np.arange(1, number_of_stages + 1)

Elements = ["Y", "Ce", "Nd", "Sm", "Gd", "Dy"]

m.pH = Var(m.fs.time, stage_number)

@m.Constraint(m.fs.time, stage_number)
def pH_value(m,t,s):
    #return m.fs.solex.mscontactor.aqueous[t,s].conc_mol_comp['H'] == 10**(-m.pH[t,s]) * units.mol/units.L
    return m.pH[t,s] == -log10(m.fs.solex.mscontactor.aqueous[t,s].conc_mol_comp['H'] * units.L/units.mol)

@m.Constraint(m.fs.time, stage_number, Elements)
def distribution_calculation(m,t,s,e):
    a, b = D_calculation(e,5)
    return m.fs.solex.distribution_coefficient[t, s, "aqueous", "organic", e] == 10**(a*m.pH[t,s] + b)


for s in stage_number:
    if s == 1:
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 5.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 3 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 24.7 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 99.1 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 32.4 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 58.2 / 100

    else:
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 4.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 12.3 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 6.4 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 16.7 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 23.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 15.1 / 100


"""
Fixation of the inlet conditions and the initial state values for all the components.

"""

m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H2O"].fix(1e6)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H"].fix(10.755)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["SO4"].fix(10)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["HSO4"].fix(1e4)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Al"].fix(422.375)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ca"].fix(109.542)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Fe"].fix(688.266)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sc"].fix(0.032)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Y"].fix(0.124)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["La"].fix(0.986)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ce"].fix(2.277)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Pr"].fix(0.303)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Nd"].fix(0.946)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sm"].fix(0.097)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Gd"].fix(0.2584)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Dy"].fix(0.047)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Cl"].fix(1e-8)

m.fs.solex.mscontactor.aqueous_inlet_state[0].flow_vol.fix(62.01)

m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Al"].fix(1.267e-5)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ca"].fix(2.684e-5)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Fe"].fix(2.873e-6)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sc"].fix(1.734)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Y"].fix(2.179e-5)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["La"].fix(0.000105)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ce"].fix(0.00031)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Pr"].fix(3.711e-5)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Nd"].fix(0.000165)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sm"].fix(1.701e-5)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Gd"].fix(3.357e-5)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Dy"].fix(8.008e-6)

m.fs.solex.mscontactor.organic_inlet_state[0].flow_vol.fix(62.01)

class AqueousPropertyScale(CustomScalerBase):

    DEFAULT_SCALING_FACTORS = {"flow_vol":1, "conc_mass_comp":1}

    def variable_scaling_routine(self, model, overwrite: bool = False, submodel_scalers: dict = None):
        self.scale_variable_by_default(model.flow_vol, overwrite=overwrite)
        for k, v in model.conc_mass_comp.items():
            if k == "H2O":
                self.set_variable_scaling_factor(v, 1, overwrite=overwrite)
            elif k == ['H','SO4','HSO4']:
                self.set_variable_scaling_factor(v, 1, overwrite=overwrite)
            else:
                self.scale_variable_by_default(v,overwrite=False)
    
    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        if model.is_property_constructed("h2o_concentration"):
            for v in model.h2o_concentration.values():
                self.scale_constraint_by_nominal_value(v, scheme="inverse_maximum", overwrite=overwrite)
        if model.is_property_constructed("molar_concentration_constraint"):
            for v in model.molar_concentration_constraint.values():
                self.scale_constraint_by_nominal_value(v, scheme="inverse_maximum", overwrite=overwrite)
        if model.is_property_constructed("hso4_dissociation"):
            for v in model.molar_concentration_constraint.values():
                self.scale_constraint_by_nominal_value(v, scheme="inverse_maximum", overwrite=overwrite)



class OrganicPropertyScale(CustomScalerBase):

    DEFAULT_SCALING_FACTORS = {"flow_vol":1, "conc_mass_comp":1}

    def variable_scaling_routine(self, model, overwrite: bool = False, submodel_scalers: dict = None):
        self.scale_variable_by_default(model.flow_vol, overwrite=overwrite)
        for k, v in model.conc_mass_comp.items():
            self.scale_variable_by_default(v,overwrite=False)
    
    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        if model.is_property_constructed("molar_concentration_constraint"):
            for v in model.molar_concentration_constraint.values():
                self.scale_constraint_by_nominal_value(v, scheme="inverse_maximum", overwrite=overwrite)


class SXScale(CustomScalerBase):

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):

        self.call_submodel_scaler_method(
            model=model,
            submodel="mscontactor.aqueous_inlet_state",
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # self.propagate_state_scaling(
        #     target_state=model.mscontactor.aqueous,
        #     source_state=model.mscontactor.aqueous_inlet_state,
        #     overwrite=overwrite,
        # )

        self.call_submodel_scaler_method(
            model=model,
            submodel="mscontactor.aqueous",
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        self.call_submodel_scaler_method(
            model=model,
            submodel="mscontactor.organic_inlet_state",
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # self.propagate_state_scaling(
        #     target_state=model.mscontactor.organic,
        #     source_state=model.mscontactor.organic_inlet_state,
        #     overwrite=overwrite,
        # )

        self.call_submodel_scaler_method(
            model=model,
            submodel="mscontactor.organic",
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        
        for v in model.mscontactor.aqueous_inherent_reaction_extent.values():
            self.set_variable_scaling_factor(v, 1e2)
        
        for v in model.mscontactor.material_transfer_term.values():
            self.set_variable_scaling_factor(v, 1e2)
    
    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):

        self.call_submodel_scaler_method(
            model=model,
            submodel="mscontactor.aqueous_inlet_state",
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # self.propagate_state_scaling(
        #     target_state=model.mscontactor.aqueous,
        #     source_state=model.mscontactor.aqueous_inlet_state,
        #     overwrite=overwrite,
        # )

        self.call_submodel_scaler_method(
            model=model,
            submodel="mscontactor.aqueous",
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        self.call_submodel_scaler_method(
            model=model,
            submodel="mscontactor.organic_inlet_state",
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # self.propagate_state_scaling(
        #     target_state=model.mscontactor.organic,
        #     source_state=model.mscontactor.organic_inlet_state,
        #     overwrite=overwrite,
        # )

        self.call_submodel_scaler_method(
            model=model,
            submodel="mscontactor.organic",
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        
        for c in model.mscontactor.component_data_objects(
            Constraint, descend_into=False
        ):
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        
        if hasattr(model, "mass_transfer_constraint"):
            for c in model.mass_transfer_constraint.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )
        
scaler = SXScale()
scaler.scale_model(
    m.fs.solex,
    submodel_scalers={
        "mscontactor.aqueous_inlet_state": AqueousPropertyScale,
        "mscontactor.aqueous": AqueousPropertyScale,
        "mscontactor.organic_inlet_state": OrganicPropertyScale,
        "mscontactor.organic": OrganicPropertyScale,
    },
)

report_scaling_factors(m, descend_into=True)

"""
# Solution of the model and display of the final results.

# """
solver = get_solver(
    "ipopt_v2", writer_config={"linear_presolve": True, "scale_model": True}
)
solver.solve(m, tee=True)

# solver = SolverFactory("ipopt")
# solver.solve(m, tee=True)

# # Final organic outlet display
# m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp.display()
# m.fs.solex.mscontactor.organic[0, 1].conc_mol_comp.display()

# # Final aqueous outlets display
# m.fs.solex.mscontactor.aqueous[0, number_of_stages].conc_mass_comp.display()
# m.fs.solex.mscontactor.aqueous[0, number_of_stages].conc_mol_comp.display()
