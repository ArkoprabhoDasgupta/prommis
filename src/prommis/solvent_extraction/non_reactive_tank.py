from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.environ import Reference

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)

__author__ = "Arkoprabho Dasgupta"


@declare_process_block_class("NonReactiveTank")
class NonReactiveTankData(UnitModelBlockData):
    """
    Standard Non reactive tank Unit Model Class
    """

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
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
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

    def build(self):
        """
        Begin building model (pre-DAE transformation).
        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(NonReactiveTankData, self).build()

        # Build Control Volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_geometry()

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_rate_reactions=False,
            has_equilibrium_reactions=False,
            has_phase_equilibrium=False,
        )

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

        # Add object references
        self.volume = Reference(self.control_volume.volume[:])
