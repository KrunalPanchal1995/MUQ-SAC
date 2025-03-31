#!/usr/bin/python
from __future__ import division
from __future__ import print_function
import pandas as pd
import numpy as np
import time
import cantera as ct
from scipy.interpolate import CubicSpline, Akima1DInterpolator
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
def find_nearest(array, value):
    array = np.asarray(np.abs(np.log(array)))
    value = np.abs(np.log(value))
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]
def find_nearest_value(xi,x,y):
    index = min(range(len(x)),key=lambda i: abs(x[i] - xi))
    return y[index]
def fast_nearest_interp(xi, x, y):
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]

gas = ct.Solution('/data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/Perturbed_Mech/mechanism_15391.yaml')

reactorTemperature = 1371.0 #Kelvin
reactorPressure = 360000.0

global criteria
if reactorTemperature> 2000:
    criteria = reactorTemperature + 20
else:
    criteria = reactorTemperature + 50

gas.TPX = reactorTemperature, reactorPressure,{'H2': 0.003, 'CO': 0.0562, 'O2': 0.0296, 'AR': 0.9112}
r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
reactorNetwork = ct.ReactorNet([r])
reactorNetwork.max_time_step = 0.0001
stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)
#index_a = stateVariableNames.index("OH*")

def ignitionDelay(df,pList, species,cond="max",specified_value ="None;None",exp_conc = "None"):
	if cond == "max":
		valid_indices = np.where((df["temperature"]>criteria))		
		peaks,_ = find_peaks(df[species].to_numpy()[valid_indices])
		time = df.index.to_numpy()
		tau = time[valid_indices[0][0]-1+peaks[0]]
	elif cond == "onset" and species != "p":
		time = df.index.to_numpy()
		conc = (df[species]).to_numpy()
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff((df[species]).to_numpy())
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "onset" and species == "p":
		time = df.index.to_numpy()
		conc = np.asarray(pList)
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff(conc)
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "dt-max" and species == "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(np.asarray(pList))
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "dt-max" and species != "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(df[species].to_numpy())
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "specific":
		if specified_value.split(";")[0] == None:
			raise Assertionerror("Input required for specified_value in ignition delay")
		else:
			target = float(specified_value.split(";")[0])
			unit = specified_value.split(";")[1]
			if unit == "molecule":
				avogadro = 6.02214E+23
				target = target/avogadro
				molecular_wt = gas.atomic_weight(species)
				target = molecular_wt*target	
				index,value = find_nearest(df[species],target)
				tau = df.index.to_numpy()[index[0]]	
			
			elif unit == "molecule/cm3":
				unit_conversion = 10**(6)
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				avogadro = 6.02214E+23
				target = (target/avogadro)*unit_conversion
				if exp_conc != "":
					exp_conc = (exp_conc/avogadro)*unit_conversion
					target = target/exp_conc
			
				f = Akima1DInterpolator(time,conc)
				time_new = np.arange(min(time),max(time),10**(-8))			
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
				
			elif unit == "mole/cm3":
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				if exp_conc != "":
					target = target/exp_conc
				
				f = Akima1DInterpolator(time,conc)
				time_new = np.arange(min(time),max(time),10**(-8))
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
	return tau

t0 = time.time()
estimatedIgnitionDelayTime = 0.04
flag = False
t = 0
pressureList = []
counter = 1;
slope_arg_max = []
count = 0
gas_T = []
gas_T.append(gas.T)
time_profile = []
time_profile.append(0)
while t<estimatedIgnitionDelayTime:
    t = reactorNetwork.step()
    gas_T.append(gas.T)
    if (counter%1 == 0):
        #if reactorNetwork.get_state().astype('float64')[index_a] <0:
         #   raise AssertionError("Numerical error!!")        	
        timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
        pressureList.append(gas.P)
    counter+=1

#while t<estimatedIgnitionDelayTime:
#    t = reactorNetwork.step()
#    #print(t,gas.T)
#    time_profile.append(t)
#    gas_T.append(gas.T)
#    time_diff = np.diff(np.asarray(time_profile))
#    conc = np.diff(np.asarray(gas_T))
#    if len(conc) == 1:
#        slope = conc/time_diff
#        index = slope[0]
#    else:
#    	slope = conc/time_diff
#    	index = int(slope.argmax())
#    slope_arg_max.append(index)
#    count+=1
#    if max(slope_arg_max) == slope_arg_max[-1] and count==1:
#    	continue
#    elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature <= 200 and t<estimatedIgnitionDelayTime:
#    	continue
#    elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature <= 200 and t>estimatedIgnitionDelayTime:
#    	break
#    elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature >= 200 and t>estimatedIgnitionDelayTime:
#    	break
#    elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature >= 100 and t<estimatedIgnitionDelayTime and flag == False:
#    	estimatedIgnitionDelayTime = t + 1e-4
#    	flag = True
#    if (counter%1 == 0):
#        #if reactorNetwork.get_state().astype('float64')[index_a] <0:
#         #   raise AssertionError("Numerical error!!")        	
#        timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
#        pressureList.append(gas.P)
#    counter+=1

tau = ignitionDelay(timeHistory,pressureList,"OH*","max","None",)
# Toc
t1 = time.time()
tau_file = open("output/tau.out",'w')
tau_file.write("#T(K)    tau(us)"+"\n"+"{}    {}".format(reactorTemperature,tau*(10**6)))
tau_file.close()
saveAll = True
if saveAll == True:
    timeHistory.to_csv("time_history.csv")
