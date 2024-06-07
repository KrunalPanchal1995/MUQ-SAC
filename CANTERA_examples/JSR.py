"""
A combustor, modeled as a single well-stirred reactor.

"""

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import os
from matplotlib import pyplot
import matplotlib.ticker as mtick


npoints = 81
T = np.linspace(700.0, 1200.0, npoints)
    
# Use reaction mechanism 
gas = ct.Solution('Syngas-NOxmech.xml') 

composition = 'H2:0.01, O2:0.01, NO:240e-6, N2:0.97976'   # phi = 0.5

states = ct.SolutionArray(gas, extra=['tres'])
residence_time = 0.24  # residence time

for i in range(npoints):
        
    gas.TPX = T[i], ct.one_atm, composition
    # Create a Reservoir for the inlet
    inlet = ct.Reservoir(gas)

    # Create the combustor, and fill it initially with a mixture consisting of the
    # equilibrium products of the inlet mixture. This state corresponds to the state
    # the reactor would reach with infinite residence time, and thus provides a good
    # initial condition from which to reach a steady-state solution on the reacting
    # branch.
    gas.equilibrate('HP')
    combustor = ct.IdealGasReactor(gas)  
    combustor.volume = 2.95e-5    # 29.5 cm3 internal volume

    # Create a reservoir for the exhaust
    exhaust = ct.Reservoir(gas)

    # Use a variable mass flow rate to keep the residence time in the reactor
    # constant (residence_time = mass / mass_flow_rate). The mass flow rate function
    # can access variables defined in the calling scope, including state variables
    # of the Reactor object (combustor) itself.

    def mdot(t):
        return combustor.mass / residence_time

    inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)

    # A PressureController has a baseline mass flow rate matching the 'master'
    # MassFlowController, with an additional pressure-dependent term. By explicitly
    # including the upstream mass flow rate, the pressure is kept constant without
    # needing to use a large value for 'K', which can introduce undesired stiffness.
    outlet_mfc = ct.PressureController(combustor, exhaust, master=inlet_mfc, K=1.0e-5)

    # the simulation only contains one reactor
    sim = ct.ReactorNet([combustor])

    sim.set_initial_time(0.0)  # reset the integrator
    sim.rtol = 1e-12
    # sim.atol = 1e-9
    sim.advance_to_steady_state(max_steps=100000)
    print('tres = {:.2e}; T = {:.1f}'.format(residence_time, combustor.T))
    states.append(combustor.thermo.state, tres=residence_time)


# Figure
plt.rc('font', family='serif')
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

X_H2, Y_H2 = [], []
for line in open('H2.txt', 'r'):
  values = [float(s) for s in line.split()]
  X_H2.append(values[0])
  Y_H2.append(values[1])

X_H2O, Y_H2O = [], []
for line in open('H2O.txt', 'r'):
  values = [float(s) for s in line.split()]
  X_H2O.append(values[0])
  Y_H2O.append(values[1])

ax1.plot(X_H2, Y_H2, linestyle='', marker='s', markersize=10, fillstyle='none', markeredgewidth=2.0, color='k')
ax1.plot(X_H2O, Y_H2O, linestyle='', marker='o', markersize=10, fillstyle='none', markeredgewidth=2.0, color='r')

ax1.legend([r'H$_2$',r'H$_2$O'], frameon=False, loc='upper right', bbox_to_anchor=(0.95, 0.95))

ax1.plot(T, states.X[:, gas.species_index('H2')], linestyle='-', linewidth=1.5, color='k')
ax1.plot(T, states.X[:, gas.species_index('H2O')], linestyle='-', linewidth=1.5, color='r')

X_NO, Y_NO = [], []
for line in open('NO.txt', 'r'):
  values = [float(s) for s in line.split()]
  X_NO.append(values[0])
  Y_NO.append(values[1])

X_NO2, Y_NO2 = [], []
for line in open('NO2.txt', 'r'):
  values = [float(s) for s in line.split()]
  X_NO2.append(values[0])
  Y_NO2.append(values[1])


ax2.plot(X_NO, Y_NO, linestyle='', marker='s', markersize=10, fillstyle='none', markeredgewidth=2.0, color='k')
ax2.plot(X_NO2, Y_NO2, linestyle='', marker='o', markersize=10, fillstyle='none', markeredgewidth=2.0, color='r')

ax2.legend([r'NO',r'NO$_2$'], frameon=False, loc='upper right', bbox_to_anchor=(0.95, 0.95))

ax2.plot(T, states.X[:, gas.species_index('NO')], linestyle='-', linewidth=1.5, color='k')
ax2.plot(T, states.X[:, gas.species_index('NO2')], linestyle='-', linewidth=1.5, color='r')

ax1.set_ylabel(r'Mole Fraction', fontsize=14)
ax1.set_xlabel(r'Temperature [K]', fontsize=14)

ax2.set_ylabel(r'Mole Fraction', fontsize=14)
ax2.set_xlabel(r'Temperature [K]', fontsize=14)
# plt.yscale('log')
ax1.set_xlim([700.0,1200.0])
ax1.set_ylim([0,0.012])

ax2.set_xlim([700.0,1200.0])
ax2.set_ylim([0,3.0e-4])
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))


width = 6.5
height = 6.0
fig1.set_size_inches(width, height)
fig2.set_size_inches(8.0, height)

fig1.suptitle(r'0.01H$_2$/0.01$_2$ in N$_2$, 240 ppm NO' '\n' r'p = 1 atm, $\phi$ = 0.5, t$_r$ = 0.24 s', fontsize=12)
fig2.suptitle(r'0.01H$_2$/0.01$_2$ in N$_2$, 240 ppm NO' '\n' r'p = 1 atm, $\phi$ = 0.5, t$_r$ = 0.24 s', fontsize=12)


fig1.savefig('phi0p5_240ppmNO_H2H2O.eps',dpi=600,format='eps')
fig2.savefig('phi0p5_240ppmNO_NONO2.eps',dpi=600,format='eps')

plt.show()