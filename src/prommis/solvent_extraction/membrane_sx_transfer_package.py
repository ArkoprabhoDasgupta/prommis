#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.common.config import ConfigValue, ConfigBlock
from pyomo.environ import Constraint, Param, Set, Var, units, log10

from idaes.core import (
    ProcessBlock,
    ProcessBlockData,
    declare_process_block_class,
    ReactionBlockBase,
    ReactionParameterBlock,
    ReactionBlockDataBase,
    MaterialFlowBasis,
)
from idaes.core.base import property_meta
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------


@declare_process_block_class("ReactionParameterTestBlock")
class _ReactionParameterBlock(ReactionParameterBlock):
    def build(self):
        super(_ReactionParameterBlock, self).build()

        self.rate_reaction_idx = Set(initialize=["r1", "r2"])
        self.equilibrium_reaction_idx = Set(initialize=["e1", "e2"])

        self.rate_reaction_stoichiometry = {
            ("r1", "p1", "c1"): 1,
            ("r1", "p1", "c2"): 1,
            ("r1", "p2", "c1"): 1,
            ("r1", "p2", "c2"): 1,
            ("r2", "p1", "c1"): 1,
            ("r2", "p1", "c2"): 1,
            ("r2", "p2", "c1"): 1,
            ("r2", "p2", "c2"): 1,
        }
        self.equilibrium_reaction_stoichiometry = {
            ("e1", "p1", "c1"): 1,
            ("e1", "p1", "c2"): 1,
            ("e1", "p2", "c1"): 1,
            ("e1", "p2", "c2"): 1,
            ("e2", "p1", "c1"): 1,
            ("e2", "p1", "c2"): 1,
            ("e2", "p2", "c1"): 1,
            ("e2", "p2", "c2"): 1,
        }

        self._reaction_block_class = ReactionBlock

        # Attribute to switch flow basis for testing
        self.basis_switch = 1

        self.set_default_scaling("reaction_rate", 101, "r1")
        self.set_default_scaling("reaction_rate", 102, "r2")

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

    @property
    def reaction_block_class(self):
        if self._reaction_block_class is not None:
            return self._reaction_block_class
        else:
            raise AttributeError(
                "{} has not assigned a ReactionBlock class to be associated "
                "with this reaction package. Please contact the developer of "
                "the reaction package.".format(self.name)
            )

    def build_reaction_block(self, *args, **kwargs):
        """
        Methods to construct a ReactionBlock associated with this
        ReactionParameterBlock. This will automatically set the parameters
        construction argument for the ReactionBlock.

        Returns:
            ReactionBlock

        """
        default = kwargs.pop("default", {})
        initialize = kwargs.pop("initialize", {})

        if initialize == {}:
            default["parameters"] = self
        else:
            for i in initialize.keys():
                initialize[i]["parameters"] = self

        return self.reaction_block_class(
            *args, **kwargs, **default, initialize=initialize
        )

    @classmethod
    def get_required_properties(cls):
        return {}


class RBlockBase(ReactionBlockBase):
    def initialize(
        blk, outlvl=idaeslog.NOTSET, optarg=None, solver=None, state_vars_fixed=False
    ):
        for k in blk.values():
            k.init_test = True


@declare_process_block_class("ReactionBlock", block_class=RBlockBase)
class ReactionBlockData(ReactionBlockDataBase):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(ReactionBlockData, self).build()

        self.reaction_rate = Var(["r1", "r2"], units=units.mol / units.m**3 / units.s)

        self.dh_rxn = {
            "r1": 10 * units.J / units.mol,
            "r2": 20 * units.J / units.mol,
            "e1": 30 * units.J / units.mol,
            "e2": 40 * units.J / units.mol,
        }

    def model_check(self):
        self.check = True

    def get_reaction_rate_basis(b):
        if b.config.parameters.basis_switch == 1:
            return MaterialFlowBasis.molar
        elif b.config.parameters.basis_switch == 2:
            return MaterialFlowBasis.mass
        else:
            return MaterialFlowBasis.other
