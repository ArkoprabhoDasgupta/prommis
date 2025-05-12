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

# Membrane_Config = ConfigDict()

# Membrane_Config.declare(
#     "property_package",
#     ConfigValue(
#         default=None,
#         domain=is_physical_parameter_block,
#         description="Property package to use for control volume",
#         doc="""Property parameter object used to define property
# calculations
# (default = 'use_parent_value')
# - 'use_parent_value' - get package from parent (default = None)
# - a ParameterBlock object""",
#     ),
# )
# Membrane_Config.declare(
#     "property_package_args",
#     ConfigValue(
#         default={},
#         description="Arguments for constructing property package",
#         doc="""A dict of arguments to be passed to the PropertyBlockData
# and used when constructing these
# (default = 'use_parent_value')
# - 'use_parent_value' - get package from parent (default = None)
# - a dict (see property package for documentation)""",
#     ),
# )
# Membrane_Config.declare(
#     "transformation_method",
#     ConfigValue(
#         default=useDefault,
#         description="Discretization method to use for DAE transformation",
#         doc="""Discretization method to use for DAE transformation. See
# Pyomo documentation for supported transformations.""",
#     ),
# )
# Membrane_Config.declare(
#     "transformation_scheme",
#     ConfigValue(
#         default=useDefault,
#         description="Discretization scheme to use for DAE transformation",
#         doc="""Discretization scheme to use when transformating domain. See
# Pyomo documentation for supported schemes.""",
#     ),
# )


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

    # CONFIG.declare(
    #     "membrane_phase",
    #     Membrane_Config(
    #         description="Membrane phase properties",
    #     ),
    # )

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

    CONFIG.declare(
        "reaction_package",
        ConfigValue(
            # TODO: Add a domain validator for this
            description="Heterogeneous reaction package for leaching.",
        ),
    )
    CONFIG.declare(
        "reaction_package_args",
        ConfigValue(
            default=None,
            domain=dict,
            description="Arguments for heterogeneous reaction package for leaching.",
        ),
    )

    def build(self):
        super().build()

        self.module_length = Var(
            domain=PositiveReals,
            initialize=1,
            doc="Bed length",
            units=units.m,
        )

        self.tube_inner_radius = Var(
            domain=PositiveReals,
            initialize=1,
            doc="Tube inner radius",
            units=units.m,
        )

        self.tube_outer_radius = Var(
            domain=PositiveReals,
            initialize=1,
            doc="Tube outer radius",
            units=units.m,
        )

        self.shell_radius = Var(
            domain=PositiveReals,
            initialize=1,
            doc="Shell radius",
            units=units.m,
        )

        self.number_of_tubes = Var(
            domain=PositiveIntegers, initialize=1, doc="Number of tubes in the module"
        )

        # Feed phase

        self.feed_phase = ControlVolume1DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.feed_phase.property_package,
            property_package_args=self.config.feed_phase.property_package_args,
            reaction_package=self.config.reaction_package,
            reaction_package_args=self.config.reaction_package_args,
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
            has_rate_reactions=True,
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
            reaction_package=self.config.reaction_package,
            reaction_package_args=self.config.reaction_package_args,
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
            has_rate_reactions=True,
        )

        self.strip_phase.add_energy_balances(
            balance_type=self.config.strip_phase.energy_balance_type,
        )

        self.strip_phase.apply_transformation()

        def shell_area(b):
            return self.strip_phase.area == pi * (
                b.shell_radius**2 - b.number_of_tubes * (b.tube_outer_radius**2)
            )

        self.shell_area_rule = Constraint(
            doc="Cross-sectional area of the tube",
            rule=shell_area,
        )

        self.add_inlet_port(name="strip_phase_inlet", block=self.strip_phase)
        self.add_outlet_port(name="strip_phase_outlet", block=self.strip_phase)

        # Membrane phase

        # self.r = ContinuousSet(
        #     bounds=(
        #         value(self.tube_outer_radius),
        #         value(self.shell_radius),
        #     )
        # )

        # # Membrane phase

        # self.membrane_phase = self.config.membrane_phase[
        #     "property_package"
        # ].build_state_block(
        #     self.flowsheet().time,
        #     self.r,
        #     self.config.control_volume_rfaces,
        #     doc="Feed phase state block",
        #     **self.config.membrane_phase["property_package_args"],
        # )
