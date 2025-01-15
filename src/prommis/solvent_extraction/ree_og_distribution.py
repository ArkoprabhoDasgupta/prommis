#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Initial property package for the organic phase solution of the solvent extraction
unit operation.

Authors: Arkoprabho Dasgupta

"""

from pyomo.environ import Param, Set, Var, units, Constraint

from idaes.core import (
    Component,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.misc import add_object_reference


@declare_process_block_class("REESolExOgParameters")
class REESolExOgParameterData(PhysicalParameterBlock):
    """
    This is a property package for the organic phase solution of the solvent extraction
    unit operation of the University of Kentucky pilot plant flowsheet.

    This  includes the following components:

    * Solvent: DEHPA
    * Rare Earths: Sc, Y, La, Ce, Pr, Nd, Sm, Gd, Dy
    * Impurities: Al, Ca, Fe

    DEHPA is not considered to be involved in any reaction.

    """

    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvents
        self.DEHPA = Component()

        # Contaminants
        self.Al = Component()
        self.Ca = Component()
        self.Fe = Component()

        # REEs
        self.Sc = Component()
        self.Y = Component()
        self.La = Component()
        self.Ce = Component()
        self.Pr = Component()
        self.Nd = Component()
        self.Sm = Component()
        self.Gd = Component()
        self.Dy = Component()

        # self.dissolved_elements = Set(
        #     initialize=[
        #         "Al",
        #         "Ca",
        #         "Fe",
        #         "Sc",
        #         "Y",
        #         "La",
        #         "Ce",
        #         "Pr",
        #         "Nd",
        #         "Sm",
        #         "Gd",
        #         "Dy",
        #     ]
        # )

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "DEHPA": 322.431e-3,
                "Sc": 44.946e-3,
                "Y": 88.905e-3,
                "La": 138.905e-3,
                "Ce": 140.116e-3,
                "Pr": 140.907e-3,
                "Nd": 144.242e-3,
                "Sm": 150.36e-3,
                "Gd": 157.25e-3,
                "Dy": 162.50e-3,
                "Al": 26.982e-3,
                "Ca": 40.078e-3,
                "Fe": 55.845e-3,
            },
        )

        # density of DEHPA
        self.dens_mass = Param(
            initialize=975.8e-3,
            units=units.kg / units.litre,
            mutable=True,
        )

        self._state_block_class = REESolExOgStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "conc_mol_comp": {"method": None},
                "dens_mass": {"method": "_dens_mass"},
            }
        )
        obj.add_default_units(
            {
                "time": units.hour,
                "mass": units.kg,
                "amount": units.mol,
                "length": units.m,
                "temperature": units.K,
            }
        )


class _REESolExOgStateBlock(StateBlock):
    def fix_initialization_states(self):
        fix_state_vars(self)


@declare_process_block_class("REESolExOgStateBlock", block_class=_REESolExOgStateBlock)
class REESolExOgStateBlockData(StateBlockData):
    """
    State block for organic phase solution of the solvent extraction process.

    """

    def build(self):
        super().build()

        self.conc_mass_comp = Var(
            self.params.component_list,
            units=units.mg / units.L,
            initialize=1e-8,
            bounds=(1e-20, None),
        )

        self.flow_vol = Var(units=units.L / units.hour, bounds=(1e-8, None))

        self.conc_mol_comp = Var(
            self.params.component_list,
            units=units.mol / units.L,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

        # Concentration conversion constraint
        @self.Constraint(self.params.component_list)
        def molar_concentration_constraint(b, j):
            return (
                units.convert(
                    b.conc_mol_comp[j] * b.params.mw[j], to_units=units.mg / units.litre
                )
                == b.conc_mass_comp[j]
            )

        if not self.config.defined_state:
            # Concentration of H2O based on assumed density
            self.dehpa_concentration = Constraint(
                expr=self.conc_mass_comp["DEHPA"] == 975.8e3 * units.mg / units.L
            )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def _dens_mass(self):
        add_object_reference(self, "dens_mass", self.params.dens_mass)

    def get_material_flow_terms(self, p, j):
        if j == "DEHPA":
            return units.convert(
                self.flow_vol * self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.hour,
            )
        else:
            return units.convert(
                self.flow_vol * self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.hour,
            )

    def get_material_density_terms(self, p, j):
        if j == "DEHPA":
            return units.convert(
                self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.m**3,
            )
        else:
            return units.convert(
                self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.m**3,
            )

    def define_state_vars(self):
        return {"flow_vol": self.flow_vol, "conc_mass_comp": self.conc_mass_comp}
