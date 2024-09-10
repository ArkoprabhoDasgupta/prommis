from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool, ConfigDict
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

@declare_process_block_class("NeutralizationTank")
class NeutralizationTankData(UnitModelBlockData):

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
        **default** = False. Equilibrium Reactors do not support dynamic behavior.""",
                ),
            )
    
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
        **default** - False. Equilibrium reactors do not have defined volume, thus
        this must be False.""",
                ),
            )
    
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

    base_config = ConfigDict()

    base_config.declare(
        "base_flowrate",
        ConfigValue(
            domain=float,
            description="Volumetric flowrate of the base"
        )
    )

    base_config.declare(
        "base_concentration",
        ConfigValue(
            domain=float,
            description="Concentration of the base"
        )
    )

    CONFIG.declare(
        "base_stream",
        ConfigDict(
            implicit=True,
            implicit_domain=base_config,
            description="Total base flowrate"
        )
    )

    def build(self):

        """
        Begin building model.

        Args:
            None

        Returns:
            None
        """

        super(NeutralizationTankData, self).build

        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        




