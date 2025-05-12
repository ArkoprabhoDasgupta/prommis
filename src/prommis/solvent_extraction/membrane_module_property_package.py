#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Initial property package for membrane based solvent extraction module.

Authors: Arkoprabho Dasgupta
"""

from pyomo.environ import (
    Constraint,
    Param,
    Set,
    Var,
    units,
    PositiveReals,
    Reals,
)

from idaes.core import (
    Component,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
    EnergyBalanceType,
)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
# Membrane solvent extraction property package
@declare_process_block_class("MembraneSXModuleParameters")
class MembraneSXModuleParameterData(PhysicalParameterBlock):

    def build(self):
        super().build()

        self.liquid = Phase()

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

        # Contaminants
        self.Al = Component()
        self.Ca = Component()
        self.Fe = Component()

        # self.mw = Param(
        #     self.component_list,
        #     units=units.kg / units.mol,
        #     initialize={
        #         "Sc": 44.946e-3,
        #         "Y": 88.905e-3,
        #         "La": 138.905e-3,
        #         "Ce": 140.116e-3,
        #         "Pr": 140.907e-3,
        #         "Nd": 144.242e-3,
        #         "Sm": 150.36e-3,
        #         "Gd": 157.25e-3,
        #         "Dy": 162.50e-3,
        #         "Al": 26.982e-3,
        #         "Ca": 40.078e-3,
        #         "Fe": 55.845e-3,
        #     },
        # )

        self.D_coeff = Param(
            self.component_list,
            units=(units.m**2) / units.s,
            initialize={
                "Sc": 10,
                "Y": 10,
                "La": 10,
                "Ce": 10,
                "Pr": 10,
                "Nd": 10,
                "Sm": 10,
                "Gd": 10,
                "Dy": 10,
                "Al": 10,
                "Ca": 10,
                "Fe": 10,
            },
        )

        self._state_block_class = MembraneSXModuleStateBlock

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
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class _MembraneSXModuleStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        # Deactivate inherent reactions
        for sbd in self.values():
            if not sbd.config.defined_state:
                sbd.h2o_concentration.deactivate()


@declare_process_block_class(
    "MembraneSXModuleStateBlock",
    block_class=_MembraneSXModuleStateBlock,
)
class MembraneSXModuleStateBlockData(StateBlockData):
    """
    State block for leach solution of West Kentucky No. 13 coal by H2SO4.

    """

    def build(self):
        super().build()

        self.membrane_flux = Var(
            self.params.component_list,
            units=units.mol / (units.sec * units.m**2),
            bounds=(1e-20, None),
        )

        self.conc_mol_comp = Var(
            self.params.component_list,
            units=units.mol / units.L,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

    def get_diffusive_transport_terms(self, p, j):

        area = self.parent_block().area
        # Note conversion to mol/hour
        if j == "H2O":
            # Assume constant density of 1 kg/L
            return units.convert(
                area * self.flow_velocity * self.params.dens_mass / self.params.mw[j],
                to_units=units.mol / units.hour,
            )

        else:
            # Need to convert from moles to mass
            return units.convert(
                area * self.flow_velocity * self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.hour,
            )

    def get_material_density_terms(self, p, j):
        if j == "H2O":
            return units.convert(
                self.params.dens_mass / self.params.mw[j],
                to_units=units.mol / units.m**3,
            )
        else:
            return units.convert(
                self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.m**3,
            )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_velocity": self.flow_velocity,
            "conc_mass_comp": self.conc_mass_comp,
        }
