#!/usr/bin/python
from __future__ import division
from __future__ import print_function
import pandas as pd
import numpy as np
import numpy
import time
import cantera as ct
import scipy.integrate
from scipy.interpolate import CubicSpline, Akima1DInterpolator
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
global dt
dt = 1e-05
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


def getTimeProfile(t1,t2):
	time =np.arange(t1,t2,1e-05)
	return time

	
def getPressureProfile(gas,time,dpdt):
	p = gas.P
	p_new = p
	PressureProfile = []
	dt = 1e-05
	for t in time:
		PressureProfile.append(p_new)
		p_new = p_new + dt*(0.01*dpdt*p_new*1000)
	return np.asarray(PressureProfile)

def getVolumeProfile_From_PressureProfile(gas,PressureProfile ):
	rho_o,Po = gas.DP
	gamma = gas.cp/gas.cv
	VolumeProfile = []
	for P_t in PressureProfile:
		VolumeProfile.append((1/rho_o)*(P_t/Po)**(-1/gamma))
	return np.asarray(VolumeProfile)
	
class VolumeProfile(object):
    def __init__(self, keywords):
        self.time = np.array(keywords["vproTime"])
        self.volume = np.array(keywords["vproVol"])/keywords["vproVol"][0]
        self.velocity = np.diff(self.volume)/np.diff(self.time)
        self.velocity = np.append(self.velocity, 0)
    def __call__(self, t):
        if t < self.time[-1]:
            prev_time_point = self.time[self.time <= t][-1]
            index = np.where(self.time == prev_time_point)[0][0]
            return self.velocity[index]
        else:
            return 0
            

###========================================================================================######
### Main code starts #####
###========================================================================================######
gas = ct.Solution('mechanism_0.yaml')

reactorTemperature = 925.0 #Kelvin
reactorPressure = 186438.0

global criteria
if reactorTemperature> 2000:
    criteria = reactorTemperature + 20
else:
    criteria = reactorTemperature + 50

gas.TPX = reactorTemperature, reactorPressure,{'H2': 0.00545, 'CO': 0.01537, 'O2': 0.01041, 'N2': 0.02876, 'AR': 0.94}
r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
env = ct.Reservoir(ct.Solution('/data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/Perturbed_Mech/mechanism_0.yaml'))

dpdt = 4.22#%/ms
estimatedIgnitionDelayTime = 0.01
flag = False
time = getTimeProfile(0,estimatedIgnitionDelayTime)
pressure_profile = getPressureProfile(gas,time,dpdt)
volume_profile = getVolumeProfile_From_PressureProfile(gas,pressure_profile)
string = "time(s),volume(cm3)\n"
for i, t in enumerate(time):
	string+=f"{t},{volume_profile[i]}\n"
g = open(f"VTIM_P_{int(reactorPressure/100000)}_T_{int(reactorTemperature)}.csv","w").write(string)

keywords = {"vproTime": time, "vproVol": volume_profile}
ct.Wall(r, env, velocity=VolumeProfile(keywords));

reactorNetwork = ct.ReactorNet([r])

reactorNetwork.max_time_step = 0.0001

stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)


####============================#####
### Solver specific information #####
###    Do not change            #####
####============================#####

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
tau = ignitionDelay(timeHistory,pressureList,"OH*","dt-max","None",)
# Toc
tau_file = open("output/tau.out",'w')
tau_file.write("#T(K)    tau(us)"+"\n"+"{}    {}".format(reactorTemperature,tau*(10**6)))
tau_file.close()
species = "OH*"
time = timeHistory.index.to_numpy()
if species != "p":
	conc = (timeHistory[species]).to_numpy()
else:
	conc = pressureList
	
fig = plt.figure()
plt.plot(time,conc,"b-",label="OH* profile")
plt.legend()
plt.savefig("profile.pdf")
saveAll = False
if saveAll == True:
    timeHistory.to_csv("time_history.csv")
