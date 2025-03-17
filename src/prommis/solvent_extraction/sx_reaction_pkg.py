#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Simple leaching model for West Kentucky No. 13 coal refuse in H2SO4.

Authors: Andrew Lee

This is an example of how to write a custom heterogeneous reaction package for use with the
LeachTrain unit model.

"""
from pyomo.common.config import ConfigValue
from pyomo.environ import Constraint, Param, Set, Var, units, log10

from idaes.core import ProcessBlock, ProcessBlockData, declare_process_block_class
from idaes.core.base import property_meta
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("SolventExtractionReactions")
class SolventExtractionReactionsData(
    ProcessBlockData, property_meta.HasPropertyClassMetadata
):

    def build(self):
        super().build()

        self._reaction_block_class = SolventExtractionReactionsBlock

        REE_list = ["La", "Y", "Pr", "Ce", "Nd", "Sm", "Gd", "Dy"]
        Impurity_list = ["Al", "Ca", "Fe", "Sc"]

        index_list = [f"{e}_mass_transfer" for e in (REE_list + Impurity_list)]
        element_list = REE_list + Impurity_list

        self.element_list = Set(initialize=element_list)
        self.reaction_idx = Set(initialize=index_list)

        reaction_stoichiometry = {}

        for e in REE_list:
            reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", e)] = -1
            reaction_stoichiometry[(f"{e}_mass_transfer", "organic", f"{e}_o")] = 1
            reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", "H")] = 3
            reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", "HSO4")] = 0
            reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", "SO4")] = 0
            reaction_stoichiometry[(f"{e}_mass_transfer", "organic", "DEHPA")] = -3

        for e in Impurity_list:
            if e == "Ca":
                reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", e)] = -1
                reaction_stoichiometry[(f"{e}_mass_transfer", "organic", f"{e}_o")] = 1
                reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", "H")] = 2
                reaction_stoichiometry[(f"{e}_mass_transfer", "organic", "DEHPA")] = -2
            else:
                reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", e)] = -1
                reaction_stoichiometry[(f"{e}_mass_transfer", "organic", f"{e}_o")] = 1
                reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", "H")] = 3
                reaction_stoichiometry[(f"{e}_mass_transfer", "organic", "DEHPA")] = -3

        self.reaction_stoichiometry = reaction_stoichiometry

        self.extractant_dosage = Param(
            doc="Extractant dosage of the system", initialize=1, mutable=True
        )

        self.m0 = Param(
            self.element_list,
            initialize={
                "Ce": 0.30916,
                "Y": 1.63166,
                "Gd": 1.0225,
                "Dy": 1.70783,
                "Sm": 0.81233,
                "Nd": 0.31183,
                "La": 0.54,
                "Pr": 0.29,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.m1 = Param(
            self.element_list,
            initialize={
                "Ce": 0.04816,
                "Y": 0.15166,
                "Gd": 0.0195,
                "Dy": 0.06443,
                "Sm": -0.02247,
                "Nd": 0.03763,
                "La": 0,
                "Pr": 0,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.B0 = Param(
            self.element_list,
            initialize={
                "Ce": -1.66021,
                "Y": -2.12601,
                "Gd": -2.24143,
                "Dy": -2.42226,
                "Sm": -2.12172,
                "Nd": -1.62372,
                "La": -1.93,
                "Pr": -1.48,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.B1 = Param(
            self.element_list,
            initialize={
                "Ce": -0.38599,
                "Y": 0.26612,
                "Gd": 0.03065,
                "Dy": -0.02538,
                "Sm": 0.17414,
                "Nd": -0.38096,
                "La": 0,
                "Pr": 0,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.K1 = Param(
            self.element_list,
            initialize={
                "Ce": 0,
                "Y": 0,
                "Gd": 0,
                "Dy": 0,
                "Sm": 0,
                "Nd": 0,
                "La": 0,
                "Pr": 0,
                "Sc": 0.991,
                "Al": 0.052,
                "Ca": 0.03,
                "Fe": 0.247,
            },
        )

        self.K2 = Param(
            self.element_list,
            initialize={
                "Ce": 0,
                "Y": 0,
                "Gd": 0,
                "Dy": 0,
                "Sm": 0,
                "Nd": 0,
                "La": 0,
                "Pr": 0,
                "Sc": 0.167,
                "Al": 0.049,
                "Ca": 0.123,
                "Fe": 0.064,
            },
        )

        self.K_corr = Param(
            self.element_list,
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

    @classmethod
    def define_metadata(cls, obj):
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

        return self.reaction_block_class(  # pylint: disable=not-callable
            *args, **kwargs, **default, initialize=initialize
        )


class _SolventExtractionReactionsBlock(ProcessBlock):
    pass


@declare_process_block_class(
    "SolventExtractionReactionsBlock", block_class=_SolventExtractionReactionsBlock
)
class SolventExtractionReactionsData(ProcessBlockData):
    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "parameters",
        ConfigValue(
            # TODO
            # domain=is_reaction_parameter_block,
            description="""A reference to an instance of the Heterogeneous Reaction Parameter
    Block associated with this property package.""",
        ),
    )

    def build(self):
        """
        Reaction block for leaching of West Kentucky No. 13 coal refuse in H2SO4.

        """
        super().build()

        add_object_reference(self, "_params", self.config.parameters)

        self.distribution_coefficient = Var(
            self.params.element_list,
            initialize=1,
        )

        def distribution_eq(b, e):
            t, s = b.index()
            aq_block = b.parent_block().aqueous[t, s]
            aq_in_block = b.parent_block().aqueous_inlet_state[t]
            org_block = b.parent_block().organic[t, s]
            org_in_block = b.parent_block().organic_inlet_state[t]

            pH = aq_block.pH_phase["liquid"]
            # dosage = b.parent_block().parent_block().extractant_dosage
            R = aq_in_block.flow_vol / org_in_block.flow_vol
            C_in = org_in_block.conc_mol_comp[f"{e}_o"] / aq_in_block.conc_mol_comp[e]

            if s == 1:
                return (b.distribution_coefficient[e]) == 10 ** (
                    (b.params.m0[e] + b.params.extractant_dosage * b.params.m1[e]) * pH
                    + (
                        b.params.B0[e]
                        + b.params.B1[e] * log10(b.params.extractant_dosage)
                    )
                ) * (1 - b.params.K_corr[e]) + b.params.K_corr[e] * (
                    (R * b.params.K1[e] + C_in) / (1 - b.params.K1[e])
                )

            else:
                return (b.distribution_coefficient[e]) == 10 ** (
                    (b.params.m0[e] + b.params.extractant_dosage * b.params.m1[e]) * pH
                    + (
                        b.params.B0[e]
                        + b.params.B1[e] * log10(b.params.extractant_dosage)
                    )
                ) * (1 - b.params.K_corr[e]) + b.params.K_corr[e] * (
                    (R * b.params.K2[e] + C_in) / (1 - b.params.K2[e])
                )

        self.distribution_constraint = Constraint(
            self.params.element_list, rule=distribution_eq
        )

    @property
    def params(self):
        return self._params
