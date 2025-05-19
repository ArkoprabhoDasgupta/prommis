from pyomo.common.config import ConfigValue, Bool, ListOf, ConfigDict, In
from pyomo.dae import DerivativeVar, ContinuousSet
from pyomo.environ import (
    Constraint,
    Param,
    Block,
    units,
    Set,
    Var,
    Reals,
    PositiveIntegers,
    PositiveReals,
    value,
    log,
)
from math import pi

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    FlowDirection,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    ControlVolume1DBlock,
)
from idaes.core.util.config import is_physical_parameter_block, DefaultBool
from idaes.core.util.constants import Constants

from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog

Stream_Config = ConfigDict()

Stream_Config.declare(
    "material_balance_type",
    ConfigValue(
        default=MaterialBalanceType.useDefault,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
    ),
)
Stream_Config.declare(
    "energy_balance_type",
    ConfigValue(
        default=EnergyBalanceType.useDefault,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
    ),
)
Stream_Config.declare(
    "momentum_balance_type",
    ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should
be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
    ),
)
Stream_Config.declare(
    "has_pressure_change",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
    ),
)
Stream_Config.declare(
    "has_phase_equilibrium",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Phase equilibrium term construction flag",
        doc="""Argument to enable phase equilibrium.
- True - include phase equilibrium term
- False - do not include phase equilibrium term""",
    ),
)
Stream_Config.declare(
    "property_package",
    ConfigValue(
        default=None,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property
calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object""",
    ),
)
Stream_Config.declare(
    "property_package_args",
    ConfigValue(
        default={},
        description="Arguments for constructing property package",
        doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)""",
    ),
)

Stream_Config.declare(
    "flow_direction",
    ConfigValue(
        default=FlowDirection.forward,
        domain=In(FlowDirection),
        doc="Direction of flow for stream",
        description="FlowDirection Enum indicating direction of "
        "flow for given stream. Default=FlowDirection.forward.",
    ),
)

Membrane_Config = ConfigDict()

Membrane_Config.declare(
    "property_package",
    ConfigValue(
        default=None,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property
calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object""",
    ),
)
Membrane_Config.declare(
    "property_package_args",
    ConfigValue(
        default={},
        description="Arguments for constructing property package",
        doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)""",
    ),
)


@declare_process_block_class("MembraneSolventExtraction")
class MembraneSolventExtractionData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "feed_phase",
        Stream_Config(
            description="Feed phase properties",
        ),
    )

    CONFIG.declare(
        "strip_phase",
        Stream_Config(
            description="Strip phase properties",
        ),
    )

    CONFIG.declare(
        "membrane_phase",
        Membrane_Config(
            description="Membrane phase properties",
        ),
    )

    CONFIG.declare(
        "tube_inner_radius",
        ConfigValue(
            default=None,
            description="Inner radius of tube in m",
            doc="User must define inner radius of tube",
        ),
    )

    CONFIG.declare(
        "tube_outer_radius",
        ConfigValue(
            default=None,
            description="Outer radius of tube in m",
            doc="User must define outer radius of tube",
        ),
    )

    CONFIG.declare(
        "shell_radius",
        ConfigValue(
            default=None,
            description="Radius of shell in m",
            doc="User must define radius of shell",
        ),
    )

    CONFIG.declare(
        "number_of_tubes",
        ConfigValue(
            default=1,
            domain=int,
            description="Number of tubes in the module",
            doc="""Number of tubes in the module.""",
        ),
    )

    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=20,
            domain=int,
            description="Number of finite elements length domain",
            doc="""Number of finite elements to use when discretizing length
            domain (default=20)""",
        ),
    )

    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
            discretizing length domain (default=3)""",
        ),
    )

    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default=useDefault,
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See
        Pyomo documentation for supported transformations.""",
        ),
    )

    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
            default=useDefault,
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transformating domain. See
        Pyomo documentation for supported schemes.""",
        ),
    )

    # CONFIG.declare(
    #     "reaction_package",
    #     ConfigValue(
    #         # TODO: Add a domain validator for this
    #         description="Heterogeneous reaction package for leaching.",
    #     ),
    # )
    # CONFIG.declare(
    #     "reaction_package_args",
    #     ConfigValue(
    #         default=None,
    #         domain=dict,
    #         description="Arguments for heterogeneous reaction package for leaching.",
    #     ),
    # )

    def build(self):
        super().build()

        self.module_length = Var(
            domain=PositiveReals,
            initialize=1,
            doc="Bed length",
            units=units.m,
        )

        self.tube_inner_radius = Param(
            initialize=self.config.tube_inner_radius,
            units=units.m,
            mutable=False,
            doc="Inner diameter of tube",
        )

        self.tube_outer_radius = Param(
            initialize=self.config.tube_outer_radius,
            units=units.m,
            mutable=False,
            doc="Inner diameter of tube",
        )

        self.shell_radius = Param(
            initialize=self.config.shell_radius,
            units=units.m,
            mutable=False,
            doc="Inner diameter of tube",
        )

        # Feed phase

        self.feed_phase = ControlVolume1DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.feed_phase.property_package,
            property_package_args=self.config.feed_phase.property_package_args,
            # reaction_package=self.config.reaction_package,
            # reaction_package_args=self.config.reaction_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )

        self.feed_phase.add_geometry(
            length_var=self.module_length,
            flow_direction=self.config.feed_phase.flow_direction,
        )

        self.feed_phase.add_state_blocks(
            information_flow=self.config.feed_phase.flow_direction,
            has_phase_equilibrium=self.config.feed_phase.has_phase_equilibrium,
        )

        self.feed_phase.add_material_balances(
            balance_type=self.config.feed_phase.material_balance_type,
            has_phase_equilibrium=self.config.feed_phase.has_phase_equilibrium,
            has_mass_transfer=True,
        )

        self.feed_phase.add_energy_balances(
            balance_type=self.config.feed_phase.energy_balance_type,
        )

        self.feed_phase.apply_transformation()

        self.add_inlet_port(name="feed_phase_inlet", block=self.feed_phase)
        self.add_outlet_port(name="feed_phase_outlet", block=self.feed_phase)

        def tube_area(b):
            return self.feed_phase.area == pi * (b.tube_inner_radius**2)

        self.tube_area_rule = Constraint(
            doc="Cross-sectional area of the tube",
            rule=tube_area,
        )

        # Strip phase

        self.strip_phase = ControlVolume1DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.strip_phase.property_package,
            property_package_args=self.config.strip_phase.property_package_args,
            # reaction_package=self.config.reaction_package,
            # reaction_package_args=self.config.reaction_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )

        self.strip_phase.add_geometry(
            length_var=self.module_length,
            flow_direction=self.config.strip_phase.flow_direction,
        )

        self.strip_phase.add_state_blocks(
            information_flow=self.config.strip_phase.flow_direction,
            has_phase_equilibrium=self.config.strip_phase.has_phase_equilibrium,
        )

        self.strip_phase.add_material_balances(
            balance_type=self.config.strip_phase.material_balance_type,
            has_phase_equilibrium=self.config.strip_phase.has_phase_equilibrium,
            has_mass_transfer=True,
        )

        self.strip_phase.add_energy_balances(
            balance_type=self.config.strip_phase.energy_balance_type,
        )

        self.strip_phase.apply_transformation()

        def shell_area(b):
            return self.strip_phase.area == pi * (
                b.shell_radius**2 - b.config.number_of_tubes * (b.tube_outer_radius**2)
            )

        self.shell_area_rule = Constraint(
            doc="Cross-sectional area of the tube",
            rule=shell_area,
        )

        self.add_inlet_port(name="strip_phase_inlet", block=self.strip_phase)
        self.add_outlet_port(name="strip_phase_outlet", block=self.strip_phase)

        # Membrane phase

        self.r = ContinuousSet(
            bounds=(
                value(self.tube_inner_radius),
                value(self.tube_outer_radius),
            )
        )

        # # Membrane phase

        self.membrane_phase = self.config.membrane_phase[
            "property_package"
        ].build_state_block(
            self.flowsheet().time,
            self.feed_phase.length_domain,
            self.r,
            doc="Feed phase state block",
            **self.config.membrane_phase["property_package_args"],
        )

        def membrane_concentration_profile(b, t, z, r, e):
            r_m = b.tube_outer_radius
            r_i = b.tube_inner_radius
            CM = b.membrane_phase[t, z, r].conc_mol_comp[e]
            CFDF = (
                b.feed_phase.properties[t, z].conc_mol_comp[e]
                * b.membrane_phase[t, z, value(r_i)].feed_distribution_coefficient[e]
            )
            CSDS = (
                b.strip_phase.properties[t, z].conc_mol_comp[e]
                * b.membrane_phase[t, z, value(r_m)].strip_distribution_coefficient[e]
            )
            return (CM - CFDF) * log(r_m / r_i) == (CSDS - CFDF) * log(r / r_i)

        self.membrane_concentration_profile_rule = Constraint(
            self.flowsheet().time,
            self.feed_phase.length_domain,
            self.r,
            self.config.membrane_phase["property_package"].component_list,
            rule=membrane_concentration_profile,
        )

        # def tube_REE_flux_formula(b, t, z, e):
        #     r_m = b.tube_outer_radius
        #     r_i = b.tube_inner_radius
        #     # if e in ["H2O", "HSO4", "H", "SO4", "Cl"]:
        #     #     return b.feed_phase.mass_transfer_term[t, z, "liquid", e] == 0
        #     # else:
        #     CFDF = (
        #         b.feed_phase.properties[t, z].conc_mol_comp[e]
        #         * b.membrane_phase[t, z, value(r_i)].feed_distribution_coefficient[e]
        #     )
        #     CSDS = (
        #         b.strip_phase.properties[t, z].conc_mol_comp[e]
        #         * b.membrane_phase[t, z, value(r_m)].strip_distribution_coefficient[e]
        #     )

        #     Diff_coeff = b.config.membrane_phase["property_package"].D_coeff[e]

        #     return b.feed_phase.mass_transfer_term[t, z, "liquid", e] == -(
        #         Diff_coeff * 2 * pi * (CFDF - CSDS) / log(r_m / r_i)
        #     )

        # self.tube_REE_flux_rule = Constraint(
        #     self.flowsheet().time,
        #     self.feed_phase.length_domain,
        #     self.config.membrane_phase["property_package"].component_list,
        #     rule=tube_REE_flux_formula,
        # )

        # def tube_non_REE_flux_formula(b, t, z, e):
        #     if e in ["H2O", "HSO4", "SO4", "Cl"]:
        #         return b.feed_phase.mass_transfer_term[t, z, "liquid", e] == 0
        #     elif e == "H":
        #         return b.feed_phase.mass_transfer_term[t, z, "liquid", e] == -3 * sum(
        #             b.feed_phase.mass_transfer_term[t, z, "liquid", i]
        #             for i in b.config.membrane_phase["property_package"].component_list
        #         )
        #     else:
        #         return Constraint.Skip

        # self.tube_non_REE_flux_rule = Constraint(
        #     self.flowsheet().time,
        #     self.feed_phase.length_domain,
        #     self.config.feed_phase["property_package"].component_list,
        #     rule=tube_non_REE_flux_formula,
        # )

        # # def shell_REE_flux_formula(b, t, z, e):
        # #     r_m = b.tube_outer_radius
        # #     r_i = b.tube_inner_radius
        # #     n = b.config.number_of_tubes
        # #     # if e in ["H2O", "HSO4", "H", "SO4", "Cl"]:
        # #     #     return b.strip_phase.mass_transfer_term[t, z, "liquid", e] == 0
        # #     # else:
        # #     CFDF = (
        # #         b.strip_phase.properties[t, z].conc_mol_comp[e]
        # #         * b.membrane_phase[t, z, value(r_i)].feed_distribution_coefficient[
        # #             e
        # #         ]
        # #     )
        # #     CSDS = (
        # #         b.strip_phase.properties[t, z].conc_mol_comp[e]
        # #         * b.membrane_phase[t, z, value(r_m)].strip_distribution_coefficient[
        # #             e
        # #         ]
        # #     )

        # #     Diff_coeff = b.config.membrane_phase["property_package"].D_coeff[e]

        # #     return b.strip_phase.mass_transfer_term[t, z, "liquid", e] == (
        # #         Diff_coeff * 2 * pi * n * (CFDF - CSDS) / log(r_m / r_i)
        # #     )

        # def shell_REE_flux_formula(b, t, z, e):
        #     n = b.config.number_of_tubes
        #     return (
        #         b.strip_phase.mass_transfer_term[t, z, "liquid", e]
        #         == -b.feed_phase.mass_transfer_term[t, z, "liquid", e] * n
        #     )

        # self.shell_flux_rule = Constraint(
        #     self.flowsheet().time,
        #     self.strip_phase.length_domain,
        #     self.config.membrane_phase["property_package"].component_list,
        #     rule=shell_REE_flux_formula,
        # )

        # def shell_non_REE_flux_formula(b, t, z, e):
        #     if e in ["H2O", "HSO4", "SO4", "Cl"]:
        #         return b.strip_phase.mass_transfer_term[t, z, "liquid", e] == 0
        #     elif e == "H":
        #         return b.strip_phase.mass_transfer_term[t, z, "liquid", e] == -3 * sum(
        #             b.strip_phase.mass_transfer_term[t, z, "liquid", i]
        #             for i in b.config.membrane_phase["property_package"].component_list
        #         )
        #     else:
        #         return Constraint.Skip

        # self.shell_non_REE_flux_rule = Constraint(
        #     self.flowsheet().time,
        #     self.feed_phase.length_domain,
        #     self.config.strip_phase["property_package"].component_list,
        #     rule=shell_non_REE_flux_formula,
        # )

        def tube_flux_formula(b, t, z, e):
            r_m = b.tube_outer_radius
            r_i = b.tube_inner_radius
            if e in ["H2O", "HSO4", "SO4", "Cl"]:
                return b.feed_phase.mass_transfer_term[t, z, "liquid", e] == 0
            elif e in b.config.membrane_phase["property_package"].component_list:
                CFDF = (
                    b.feed_phase.properties[t, z].conc_mol_comp[e]
                    * b.membrane_phase[t, z, value(r_i)].feed_distribution_coefficient[
                        e
                    ]
                )
                CSDS = (
                    b.strip_phase.properties[t, z].conc_mol_comp[e]
                    * b.membrane_phase[t, z, value(r_m)].strip_distribution_coefficient[
                        e
                    ]
                )

                Diff_coeff = b.config.membrane_phase["property_package"].D_coeff[e]

                return b.feed_phase.mass_transfer_term[t, z, "liquid", e] == -(
                    Diff_coeff * 2 * pi * (CFDF - CSDS) / log(r_m / r_i)
                )
            else:
                return b.feed_phase.mass_transfer_term[t, z, "liquid", e] == -3 * sum(
                    b.feed_phase.mass_transfer_term[t, z, "liquid", i]
                    for i in b.config.membrane_phase["property_package"].component_list
                )

        self.tube_flux_rule = Constraint(
            self.flowsheet().time,
            self.feed_phase.length_domain,
            self.config.feed_phase["property_package"].component_list,
            rule=tube_flux_formula,
        )

        # def tube_non_REE_flux_formula(b, t, z, e):
        #     if e in ["H2O", "HSO4", "SO4", "Cl"]:
        #         return b.feed_phase.mass_transfer_term[t, z, "liquid", e] == 0
        #     elif e == "H":
        #         return b.feed_phase.mass_transfer_term[t, z, "liquid", e] == -3 * sum(
        #             b.feed_phase.mass_transfer_term[t, z, "liquid", i]
        #             for i in b.config.membrane_phase["property_package"].component_list
        #         )
        #     else:
        #         return Constraint.Skip

        # self.tube_non_REE_flux_rule = Constraint(
        #     self.flowsheet().time,
        #     self.feed_phase.length_domain,
        #     self.config.feed_phase["property_package"].component_list,
        #     rule=tube_non_REE_flux_formula,
        # )

        # def shell_REE_flux_formula(b, t, z, e):
        #     r_m = b.tube_outer_radius
        #     r_i = b.tube_inner_radius
        #     n = b.config.number_of_tubes
        #     # if e in ["H2O", "HSO4", "H", "SO4", "Cl"]:
        #     #     return b.strip_phase.mass_transfer_term[t, z, "liquid", e] == 0
        #     # else:
        #     CFDF = (
        #         b.strip_phase.properties[t, z].conc_mol_comp[e]
        #         * b.membrane_phase[t, z, value(r_i)].feed_distribution_coefficient[
        #             e
        #         ]
        #     )
        #     CSDS = (
        #         b.strip_phase.properties[t, z].conc_mol_comp[e]
        #         * b.membrane_phase[t, z, value(r_m)].strip_distribution_coefficient[
        #             e
        #         ]
        #     )

        #     Diff_coeff = b.config.membrane_phase["property_package"].D_coeff[e]

        #     return b.strip_phase.mass_transfer_term[t, z, "liquid", e] == (
        #         Diff_coeff * 2 * pi * n * (CFDF - CSDS) / log(r_m / r_i)
        #     )

        def shell_flux_formula(b, t, z, e):
            n = b.config.number_of_tubes
            return (
                b.strip_phase.mass_transfer_term[t, z, "liquid", e]
                == -b.feed_phase.mass_transfer_term[t, z, "liquid", e] * n
            )

        self.shell_flux_rule = Constraint(
            self.flowsheet().time,
            self.strip_phase.length_domain,
            self.config.strip_phase["property_package"].component_list,
            rule=shell_flux_formula,
        )
