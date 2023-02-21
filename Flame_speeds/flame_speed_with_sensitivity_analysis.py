#!/usr/bin/env python
# coding: utf-8

# # Flame Speed with Sensitivity Analysis

# In this example we simulate a freely-propagating, adiabatic, 1-D flame and
# * Calculate its laminar burning velocity
# * Perform a sensitivity analysis of its kinetics
# 
# The figure below illustrates the setup, in a flame-fixed co-ordinate system. The reactants enter with density $\rho_{u}$, temperature $T_{u}$ and speed $S_{u}$. The products exit the flame at speed $S_{b}$, density $\rho_{b}$ and temperature $T_{b}$.

# <img src="images/flameSpeed.png" alt="Freely Propagating Flame" style="width: 300px;"/>

# ### Import Modules

# In[1]:
from __future__ import print_function
from __future__ import division
import cantera as ct
import numpy as np

print("Running Cantera Version: " + str(ct.__version__))

# In[2]:


# Import plotting modules and define plotting preference
#get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pylab as plt
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.figsize'] = (8,6)

# Get the best of both ggplot and seaborn
plt.style.use('ggplot')
plt.style.use('seaborn-deep')
plt.rcParams['figure.autolayout'] = True
# Import Pandas for DataFrames
import pandas as pd
# ### Define the reactant conditions, gas mixture and kinetic mechanism associated with the gas
# In[3]:
#Inlet Temperature in Kelvin and Inlet Pressure in Pascals
#In this case we are setting the inlet T and P to room temperature conditions
To = 914
#Po = 101325
Po  = 101325*0.047
#Define the gas-mixutre and kinetics
#In this case, we are choosing a GRI3.0 gas
#gas = ct.Solution('Davis_H2CO_2005.cti')

gas = ct.Solution('gri30.cti')
# Create a stoichiometric CH4/Air premixed mixture 
#gas.set_equivalence_ratio(3.05, {'H2':0.281,'CO':0.281},{'O2':0.092, 'N2':0.346})
#gas.X = {'H2':0.092,'CO':0.078,'O2':0.174, 'N2':0.656}
gas.X = {'H2':0.397,'CO':0.00,'O2':0.103, 'AR':0.5}
#gas.set_equivalence_ratio(1.0, 'CH4', {'O2':1.0, 'N2':3.76})
gas.TP = To, Po#,{'H2':3.05434783,'CO':3.05434783,'O2':1, 'N2':3.76}
print(gas())
# ### Define flame simulation conditions
# In[4]:
# Domain width in metres
width = 0.014
# Create the flame object
flame = ct.FreeFlame(gas, width=width)
# Define tolerances for the solver
flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
# Define logging level
loglevel = 1
# ### Solve
# In[5]:
flame.solve(loglevel=loglevel, auto=True)
Su0 = flame.velocity[0]
print("Flame Speed is: {:.2f} cm/s".format(Su0*100))
# Note that the variable Su0 will also be used downsteam in the sensitivity analysis
# ### Plot figures
# 
# Check and see if all has gone well. Plot temperature and species fractions to see
# 
# #### Temperature Plot
# In[6]:
plt.figure()
plt.plot(flame.grid*100, flame.T, '-o')
plt.xlabel('Distance (cm)')
plt.ylabel('Temperature (K)');
# #### Major species' plot
# To plot species, we first have to identify the index of the species in the array
# For this, cut & paste the following lines and run in a new cell to get the index
# 
#     for i, specie in enumerate(gas.species()):
#         print(str(i) + '. ' + str(specie))
# In[7]:
# Extract concentration data
X_CH4 = flame.X[13]
X_CO2 = flame.X[15]
X_H2O = flame.X[5]
plt.figure()
plt.plot(flame.grid*100, X_CH4, '-o', label=r'$CH_{4}$')
plt.plot(flame.grid*100, X_CO2, '-s', label=r'$CO_{2}$')
plt.plot(flame.grid*100, X_H2O, '-<', label=r'$H_{2}O$')
plt.legend(loc=2)
plt.xlabel('Distance (cm)')
plt.ylabel('MoleFractions');
plt.show()
# ## Sensitivity Analysis
# 
# See which reactions effect the flame speed the most
# In[8]:
# Create a dataframe to store sensitivity-analysis data
sensitivities = pd.DataFrame(data=[], index=gas.reaction_equations(range(gas.n_reactions)))
# ### Compute sensitivities
# In[9]:
# Set the value of the perturbation
dk = 1e-2
# Create an empty column to store the sensitivities data
sensitivities["baseCase"] = ""
# In[10]:
for m in range(gas.n_reactions):
    gas.set_multiplier(1.0) # reset all multipliers                                            
    gas.set_multiplier(1+dk, m) # perturb reaction m   
    # Always force loglevel=0 for this
    # Make sure the grid is not refined, otherwise it won't strictly 
    # be a small perturbation analysis
    flame.solve(loglevel=0, refine_grid=False)
    # The new flame speed
    Su = flame.velocity[0]
    sensitivities["baseCase"][m] = (Su-Su0)/(Su0*dk)
# This step is essential, otherwise the mechanism will have been altered
gas.set_multiplier(1.0)
# In[11]:
sensitivities.head()
# ### Make plots
# In[12]:
# Reaction mechanisms can contains thousands of elementary steps. Choose a threshold
# to see only the top few
threshold = 0.03
firstColumn = sensitivities.columns[0]
# For plotting, collect only those steps that are above the threshold
# Otherwise, the y-axis gets crowded and illegible
sensitivitiesSubset = sensitivities[sensitivities[firstColumn].abs() > threshold]
indicesMeetingThreshold = sensitivitiesSubset[firstColumn].abs().sort_values(ascending=False).index
sensitivitiesSubset.loc[indicesMeetingThreshold].plot.barh(title="Sensitivities for GRI 3.0",
                                                          legend=None)
plt.gca().invert_yaxis()
plt.rcParams.update({'axes.labelsize': 20})
plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$');
plt.show()
# Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
# plt.savefig('sensitivityPlot', dpi=300)

