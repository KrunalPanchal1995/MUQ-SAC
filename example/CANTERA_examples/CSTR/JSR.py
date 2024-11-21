import pandas as pd
import time
import cantera as ct

gas = ct.Solution("galway.yaml")
reactor_temperature = 925  # Kelvin
reactor_pressure = 1.046138 * ct.one_atm  # in atm. This equals 1.06 bars
inlet_concentrations = {"NC7H16": 0.005, "O2": 0.0275, "HE": 0.9675}
gas.TPX = reactor_temperature, reactor_pressure, inlet_concentrations
species = "CO"
# Reactor parameters
residence_time = 2  # s
reactor_volume = 30.5 * (1e-2) ** 3  # m3
max_simulation_time = residence_time+1
fuel_air_mixture_tank = ct.Reservoir(gas)
exhaust = ct.Reservoir(gas)

stirred_reactor = ct.IdealGasReactor(gas, energy="off", volume=reactor_volume)

mass_flow_controller = ct.MassFlowController(
    upstream=fuel_air_mixture_tank,
    downstream=stirred_reactor,
    mdot=stirred_reactor.mass / residence_time,
)

pressure_regulator = ct.PressureController(
    upstream=stirred_reactor, downstream=exhaust, master=mass_flow_controller
)

reactor_network = ct.ReactorNet([stirred_reactor])

time_history = ct.SolutionArray(gas, extra=["t"])

tic = time.time()

t = 0
counter = 1
while t < max_simulation_time:
    t = reactor_network.step()

    # We will store only every 10th value. Remember, we have 1200+ species, so there will be
    # 1200+ columns for us to work with
    if counter % 10 == 0:
        # Extract the state of the reactor
        time_history.append(stirred_reactor.thermo.state, t=t)

    counter += 1

# Stop the stopwatch
toc = time.time()
MF = time_history(species).X[-1][0]
tau_file = open("output/jsr.out",'w')
tau_file.write("#T(K)    mole fraction"+"\n"+"{}    {}".format(reactor_temperature,MF))
tau_file.close()  
