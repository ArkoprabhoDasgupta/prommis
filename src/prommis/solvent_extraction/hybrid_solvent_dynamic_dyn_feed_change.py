#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

""" 
Demonstration flowsheet for dynamic state solvent extraction loading process
using parameters and data derived from West Kentucky No. 13 coal refuse.

Authors: Arkoprabho Dasgupta

"""

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    units,
    TransformationFactory,
    Var,
    value,
)
from pyomo.dae.flatten import flatten_dae_components

import numpy as np

import matplotlib.pyplot as plt

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.util import from_json

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.hybrid_solvent_extraction import SolventExtraction
from prommis.solvent_extraction.D_model_systems import D_calculation

"""
Method of building a dynamic solvent extraction model with a specified number of
stages and with two separate property packages for the two inlet streams.
This is a loading operation, so no additional argument has to be specified.

"""

m = ConcreteModel()

time_duration = 60

m.fs = FlowsheetBlock(dynamic=True, time_set=[0, time_duration], time_units=units.hour)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()

number_of_stages = 3

m.fs.solex = SolventExtraction(
    number_of_finite_elements=number_of_stages,
    aqueous_stream={
        "property_package": m.fs.leach_soln,
        "flow_direction": FlowDirection.forward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
    organic_stream={
        "property_package": m.fs.prop_o,
        "flow_direction": FlowDirection.backward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
)


"""
Discretization of the time domain, and specification of the partition coefficients,
volume, volume fractions, and the initial conditions of state variables for the components
for all the stages.

"""

m.discretizer = TransformationFactory("dae.collocation")
m.discretizer.apply_to(m, nfe=12, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")

"""
Specifications of the partition coefficients, volume and volume fractions for all
the stages.

"""


m.fs.solex.mscontactor.volume[:].fix(0.4)

m.fs.solex.mscontactor.volume_frac_stream[:, :, "organic"].fix(0.4)

number_of_stages = 3
stage_number = np.arange(1, number_of_stages + 1)

Elements = ["Y", "Ce", "Nd", "Sm", "Gd", "Dy"]

"""
Initialization of the model, which gives a good starting point.

"""

from_json(m, fname="hybrid_solvent_extraction.json")


def copy_first_steady_state(m):
    # Function that propogates initial steady state guess to future time points
    # regular_vars
    regular_vars, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
    # Copy initial conditions forward
    for var in time_vars:
        for t in m.fs.time:
            if t == m.fs.time.first():
                continue
            else:
                var[t].value = var[m.fs.time.first()].value
                # var.pprint()


copy_first_steady_state(m)

# pH_loading = 1.55
m.pH_loading = Var(m.fs.time)


@m.Constraint(m.fs.time)
def pH_variance(m, t):
    if t <= 36:
        return m.pH_loading[t] == 1.3
    else:
        # return m.pH_loading[t] == 1.55 + (1.7-1.55)*(t-10)/14
        return m.pH_loading[t] == 1.4


# for e in Elements:
#     for t in m.fs.time:
#         m.fs.solex.distribution_coefficient[t,:, "aqueous", "organic", e] = D_calculation(
#             e, "5% dehpa 10% tbp", m.pH_loading[t]
#         )m.p

for e in Elements:
    for t in m.fs.time:
        if t <= 36:
            pH_loading = 1.3
            m.fs.solex.distribution_coefficient[t, :, "aqueous", "organic", e] = (
                D_calculation(e, "5% dehpa 10% tbp", pH_loading)
            )
        else:
            # pH_loading = 1.2 + (1.7-1.2)*(t-50)/10
            pH_loading = 1.41
            m.fs.solex.distribution_coefficient[t, :, "aqueous", "organic", e] = (
                D_calculation(e, "5% dehpa 10% tbp", pH_loading)
            )

for s in stage_number:
    if s == 1:
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 5.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 3 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 24.7 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 99.1 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 32.4 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 58.2 / 100

    else:
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 4.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 12.3 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 6.4 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 16.7 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 23.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 15.1 / 100


"""
Fixation of the inlet conditions and the initial state values for all the components.

"""

m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H2O"].fix(1e6)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H"].fix(1.755)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["SO4"].fix(3999.818)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["HSO4"].fix(693.459)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Al"].fix(422.375)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ca"].fix(109.542)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Cl"].fix(1e-7)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Fe"].fix(688.266)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sc"].fix(0.032)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Y"].fix(0.124)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["La"].fix(0.986)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ce"].fix(2.277)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Pr"].fix(0.303)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Nd"].fix(0.946)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sm"].fix(0.097)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Gd"].fix(0.2584)
m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Dy"].fix(0.047)

for t in m.fs.time:
    if t <= 24:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].conc_mass_comp["Y"].fix(0.124)
    else:
        m.fs.solex.mscontactor.aqueous_inlet_state[t].conc_mass_comp["Y"].fix(
            0.124 + 0.1
        )

m.fs.solex.mscontactor.aqueous_inlet_state[:].flow_vol.fix(62.01)

m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Al"].fix(1.267e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ca"].fix(2.684e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Fe"].fix(2.873e-6)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sc"].fix(1.734)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Y"].fix(2.179e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["La"].fix(0.000105)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ce"].fix(0.00031)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Pr"].fix(3.711e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Nd"].fix(0.000165)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sm"].fix(1.701e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Gd"].fix(3.357e-5)
m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Dy"].fix(8.008e-6)

m.fs.solex.mscontactor.organic_inlet_state[:].flow_vol.fix(62.01)

m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["HSO4"].fix(693.459)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["SO4"].fix(3999.818)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Al"].fix(422.375)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ca"].fix(109.542)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Cl"].fix(1e-7)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Fe"].fix(688.266)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sc"].fix(0.032)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Y"].fix(0.124)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["La"].fix(0.986)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ce"].fix(2.277)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Pr"].fix(0.303)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Nd"].fix(0.946)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sm"].fix(0.097)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Gd"].fix(0.2584)
m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Dy"].fix(0.047)

m.fs.solex.mscontactor.aqueous_inherent_reaction_extent[0, :, "Ka2"].fix(1e-8)
m.fs.solex.mscontactor.aqueous[0, :].flow_vol.fix(62.01)


m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Al"].fix(1.267e-5)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ca"].fix(2.684e-5)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Fe"].fix(2.837e-6)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sc"].fix(1.734)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["La"].fix(0.000105)
m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Pr"].fix(3.711e-5)

m.fs.solex.mscontactor.organic[0, :].flow_vol.fix(62.01)

for e in Elements:
    m.fs.solex.mscontactor.material_transfer_term[0.0, :, "aqueous", "organic", e].fix(
        1e-8
    )


"""
Solution of the model and display of the final results.

"""

solver = SolverFactory("ipopt")
solver.solve(m, tee=True)

# Final organic outlet display
m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp.display()
m.fs.solex.mscontactor.organic[time_duration, 1].conc_mol_comp.display()

# Final aqueous outlets display
m.fs.solex.mscontactor.aqueous[time_duration, number_of_stages].conc_mass_comp.display()
m.fs.solex.mscontactor.aqueous[time_duration, number_of_stages].conc_mol_comp.display()

percent_recovery = {}
for ei, e in enumerate(Elements):
    for si, s in enumerate(stage_number):
        percent_recovery[e, s] = [
            (
                1
                - (
                    (
                        m.fs.solex.mscontactor.aqueous[t, s].conc_mass_comp[e]()
                        * m.fs.solex.mscontactor.aqueous[t, s].flow_vol()
                    )
                    / (
                        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp[
                            e
                        ]()
                        * m.fs.solex.mscontactor.aqueous_inlet_state[0].flow_vol()
                    )
                )
            )
            * 100
            for t in m.fs.time
        ]


for e in Elements:
    plt.plot(m.fs.time, percent_recovery[e, 1])
plt.legend(Elements)
plt.xlabel("time, hrs")
plt.ylabel("percent recovery, %")
plt.title("Aqueous phase percent recovery graph w.r.t. time")

print(percent_recovery["Y", 1])
