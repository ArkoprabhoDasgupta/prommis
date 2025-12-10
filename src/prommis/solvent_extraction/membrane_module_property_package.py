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
    log10,
)
from pyomo.dae import DerivativeVar

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

        self.D_coeff = Param(
            self.component_list,
            units=(units.m**2) / units.s,
            initialize={
                "Sc": 1e-9,
                "Y": 1e-9,
                "La": 1e-9,
                "Ce": 1e-9,
                "Pr": 1e-9,
                "Nd": 1e-9,
                "Sm": 1e-9,
                "Gd": 1e-9,
                "Dy": 1e-9,
                "Al": 1e-9,
                "Ca": 1e-9,
                "Fe": 1e-9,
            },
        )

        self.extractant_dosage = Param(
            doc="Extractant dosage of the system", initialize=1, mutable=True
        )

        self.m0 = Param(
            self.component_list,
            initialize={
                "Ce": 0.409,
                "Y": 1.923,
                "Gd": 1.096,
                "Dy": 1.689,
                "Sm": 0.764,
                "Nd": 0.319,
                "La": 0.541,
                "Pr": 0.553,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.m1 = Param(
            self.component_list,
            initialize={
                "Ce": 0.0218,
                "Y": 0.0748,
                "Gd": 5.8e-9,
                "Dy": 0.05774,
                "Sm": 0.0218,
                "Nd": 0.00294,
                "La": 0,
                "Pr": 0,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.B0 = Param(
            self.component_list,
            initialize={
                "Ce": -1.908,
                "Y": -2.186,
                "Gd": -2.316,
                "Dy": -2.354,
                "Sm": -2.153,
                "Nd": -1.82,
                "La": -2.51,
                "Pr": -3.015,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.B1 = Param(
            self.component_list,
            initialize={
                "Ce": 0.04,
                "Y": 6.28e-7,
                "Gd": 0.208,
                "Dy": 1.09e-6,
                "Sm": 0.0086,
                "Nd": 4.32e-6,
                "La": 0.814,
                "Pr": 1.533,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.K1 = Param(
            self.component_list,
            initialize={
                "Ce": 0,
                "Y": 0,
                "Gd": 0,
                "Dy": 0,
                "Sm": 0,
                "Nd": 0,
                "La": 0,
                "Pr": 0,
                "Sc": 632.4976,
                "Al": 0.0531,
                "Ca": 0.0658,
                "Fe": 0.1496,
            },
        )

        self.K_corr = Param(
            self.component_list,
            initialize={
                "Ce": 0,
                "Y": 0,
                "Gd": 0,
                "Dy": 0,
                "Sm": 0,
                "Nd": 0,
                "La": 0,
                "Pr": 0,
                "Sc": 1,
                "Al": 1,
                "Ca": 1,
                "Fe": 1,
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
        # for sbd in self.values():
        #     if not sbd.config.defined_state:
        #         sbd.h2o_concentration.deactivate()


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

        self.feed_distribution_coefficient = Var(
            self.params.component_list,
            initialize=1,
        )

        self.strip_distribution_coefficient = Var(
            self.params.component_list,
            initialize=1,
        )

        def feed_distribution_eq(b, e):
            feed_block = b.parent_block().feed_phase.properties[
                b.index()[0], b.index()[1]
            ]

            pH = feed_block.pH_phase
            return (b.feed_distribution_coefficient[e]) == 10 ** (
                (b.params.m0[e] + b.params.extractant_dosage * b.params.m1[e]) * pH
                + (b.params.B0[e] + b.params.B1[e] * log10(b.params.extractant_dosage))
            ) * (1 - b.params.K_corr[e]) + b.params.K_corr[e] * b.params.K1[e]

        self.feed_distribution_constraint = Constraint(
            self.params.component_list, rule=feed_distribution_eq
        )

        def strip_distribution_eq(b, e):
            strip_block = b.parent_block().strip_phase.properties[
                b.index()[0], b.index()[1]
            ]

            pH = strip_block.pH_phase
            return (b.strip_distribution_coefficient[e]) == 10 ** (
                (b.params.m0[e] + b.params.extractant_dosage * b.params.m1[e]) * pH
                + (b.params.B0[e] + b.params.B1[e] * log10(b.params.extractant_dosage))
            ) * (1 - b.params.K_corr[e]) + b.params.K_corr[e] * b.params.K1[e]

        self.strip_distribution_constraint = Constraint(
            self.params.component_list, rule=strip_distribution_eq
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar
