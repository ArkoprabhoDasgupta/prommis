from pyomo.common.config import ConfigBlock, ConfigValue, In, ConfigDict
from pyomo.environ import Param, Constraint, Var, units

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
    FlowDirection,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
)


@declare_process_block_class("NeutralizationTank")
class NeutralizationTankData(UnitModelBlockData):

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
        Begin building model.

        Args:
            None

        Returns:
            None
        """

        super(NeutralizationTankData, self).build()

        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(
            information_flow=FlowDirection.forward, has_phase_equilibrium=True
        )

        self.control_volume.add_total_component_balances(
            has_rate_reactions=False,
            has_equilibrium_reactions=False,
            has_phase_equilibrium=False,
            has_mass_transfer=True,
            custom_molar_term=None,
            custom_mass_term=None,
        )

        self.base_concentration = Var(self.flowsheet().time, units=units.mol / units.L)

        self.base_flowrate = Var(self.flowsheet().time, units=units.L / units.hr)

        self.water_feed_conc = Var(self.flowsheet().time, units=units.mol / units.L)

        def water_inlet_conc(self, t):

            water_concentration = 55.55 * units.mol / units.L
            specific_volume_base = 18.8 * units.mL / units.mol

            return (
                self.water_feed_conc[t]
                == (
                    1 - units.convert(self.base_concentration[t] * specific_volume_base)
                )
                * water_concentration
            )

        self.water_feed_constraint = Constraint(
            self.flowsheet().time, rule=water_inlet_conc
        )

        def mass_transfer_term(self, t, p, c):

            if c == "H2O":
                return (
                    self.control_volume.mass_transfer_term[t, p, c]
                    == (self.base_concentration[t] + self.water_feed_conc[t])
                    * self.base_flowrate[t]
                )
            elif c == "H":
                return (
                    self.control_volume.mass_transfer_term[t, p, c]
                    == -self.base_concentration[t] * self.base_flowrate[t]
                )
            else:
                return (
                    self.control_volume.mass_transfer_term[t, p, c]
                    == 0 * units.mol / units.L
                )

        pc_set = self.config.property_package._phase_component_set

        self.mass_transfer_constraint = Constraint(
            self.flowsheet().time, pc_set, rule=mass_transfer_term
        )

        self.add_inlet_port()
        self.add_outlet_port()
