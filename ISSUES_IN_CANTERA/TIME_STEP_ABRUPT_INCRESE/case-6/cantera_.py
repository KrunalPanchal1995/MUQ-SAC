#!/usr/bin/python
from __future__ import division
from __future__ import print_function
import pandas as pd
import numpy as np
import time
import cantera as ct
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

def find_nearest(array, value):
    array = np.asarray(np.abs(np.log(array)))
    value = np.abs(np.log(value))
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

def fast_nearest_interp(xi, x, y):
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]

gas = ct.Solution('mechanism.yaml')

reactorTemperature = 1944.0 #Kelvin
reactorPressure = 112000.00000000001

criteria = 1800

if (reactorTemperature >= criteria):
	criteria = reactorTemperature + 300
elif (reactorTemperature - criteria)<=200:
	criteria = reactorTemperature + 300

gas.TPX = reactorTemperature, reactorPressure,{'H2': 0.0017, 'CO': 0.033, 'O2': 0.0347, 'Ar': 0.9306}
r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
reactorNetwork = ct.ReactorNet([r])
reactorNetwork.max_time_step = 1e-6
stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)
#index_a = stateVariableNames.index("OH*")

def ignitionDelay(df,pList, species,cond="max",specified_value ="None;None",exp_conc = "None"):
	if cond == "max":
		tau = df[species].idxmax()
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
					conc = conc*exp_conc
			
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))			
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
				
			elif unit == "mole/cm3":
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				if exp_conc != "":
					conc = conc*exp_conc
				
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
	return tau

t0 = time.time()
estimatedIgnitionDelayTime = 0.02
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
    #print(t,gas.T)
    time_profile.append(t)
    gas_T.append(gas.T)
    time_diff = np.diff(np.asarray(time_profile))
    conc = np.diff(np.asarray(gas_T))
    if len(conc) == 1:
        slope = conc/time_diff
        index = slope[0]
    else:
    	slope = conc/time_diff
    	index = int(slope.argmax())
    slope_arg_max.append(index)
    count+=1
    if max(slope_arg_max) == slope_arg_max[-1] and count==1:
    	continue
    elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature <= 200 and t<estimatedIgnitionDelayTime:
    	continue
    elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature <= 200 and t>estimatedIgnitionDelayTime:
    	break
    elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature >= 200 and t>estimatedIgnitionDelayTime:
    	break
    elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature >= 300 and t<estimatedIgnitionDelayTime and flag == False:
    	estimatedIgnitionDelayTime = 5*t
    	flag = True
    	print(flag)
    if (counter%1 == 0):
        #if reactorNetwork.get_state().astype('float64')[index_a] <0:
         #   raise AssertionError("Numerical error!!")        	
        timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
        pressureList.append(gas.P)
    counter+=1

tau = ignitionDelay(timeHistory,pressureList,"OH*","max","None",)
# Toc
t1 = time.time()
tau_file = open("output/tau.out",'w')
tau_file.write("#T(K)    tau(us)"+"\n"+"{}    {}".format(reactorTemperature,tau*(10**6)))
tau_file.close()
saveAll = True
if saveAll == True:
    timeHistory.to_csv("time_history.csv")
