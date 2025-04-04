import subprocess
import os
import yaml
from yaml import load, dump
try:
	from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
	from yaml import Loader, Dumper

def create_JPDAP_input(input_dict):
	instring = f'''* Uncertain Arrhenius parameters [A/An/AE/AnE]
{input_dict["uncertain_parameters"]}
* Uncertainty type [3slog10k/2slog10k/1slog10k/1slnk/2slnk/3slnk]
{input_dict["uncertainty_type"]}
* Number of data (1st row),data in rows: temperature uncertainty
{input_dict["len_temp_data"]}
{input_dict["temperature_unsrt_data"]}
* Test sets: number of sets (1st row),sets in rows [sa,sn,se,ran,rae,rne, omit missing] a=lnA,e=E/R
{input_dict["L"]}'''
	return instring
	
def create_SAMAP_input(input_dict):
	instring = f'''* Uncertain Arrhenius parameters [A/An/AE/AnE]
{input_dict["uncertain_parameters"]}
* Uncertainty type [3slog10k/2slog10k/1slog10k/1slnk/2slnk/3slnk]
{input_dict["uncertainty_type"]}
* Mean values of the uncertain Arrhenius parameters (A,n,E/R[K]) (skip missing parameters)
{input_dict["alpha"]} {input_dict["n"]} {input_dict["epsilon"]}
* Covariance matrix [(a,n,e)x(a,n,e)], where a=lnA, e=E/R (skip rows/columns for missing ones) 
 {input_dict["covariance_matrix"]}
* if n is uncertain limits can be set for n: n_min n_max (otherwise the next line is skipped)
{input_dict["n_min"]} {input_dict["n_max"]}
* Temperature (T) range of validity: Tmin Tmax (used for checking the uncertainty limits)
{input_dict["T_begin"]} {input_dict["T_end"]}
* Number of equidistant T points within the T range used for the discretization (>=10)
{input_dict["equidistant_T"]}
* Distribution of transformed Arrhenius parameters (UNIFORM/NORMAL)
{input_dict["sampling_distribution"]}
* Sampling method (EQUIDISTANT/RANDOM/LATIN_HC/ORTHOGONAL/SOBOL)
{input_dict["sampling_method"]}
* Random seed (integer between 1 and 2^32-1) for random number generator
{input_dict["Random_seed"]}
* Number of samples to generate and number of samples to be skipped for SOBOL (optional)
{input_dict["samples"]} {input_dict["samples_skipped"]}'''

	return instring

def create_start_profile_input(inputs,target):	
	instring = '''#########
# Input #
#########

Inputs:
 
 bin: %(bin)s
 
 Flame: %(flame_config)s
 
 MechanismFile: %(pre_file)s
 
 StartProfilesFile: %(init_profile)s
 
 ComputeWithRadiation: %(cwr)s
 
 Thermodiffusion: %(td)s

 ExpTempFile: %(etf)s
 
 CopyTo: %(copy)s
 
 UniqueID: %(id)s
 
#######################
# Boundary conditions #
#######################

Boundary:
 
 fuel: %(fuel)s
 
 oxidizer: %(oxidizer)s
 
 bathGas: %(bathGas)s
 
 globalReaction: %(glob_rxn)s
 
 pressure: %(to_Pa)s 
 
 temperature: %(to_K)s
 
 flowRate: %(flow)s
 
 units: %(units)s'''% {"bin":inputs["Bin"]["bin"],"flame_config":inputs["StartProfilesData"][target.target]["Flame"],"pre_file":inputs["Locations"]["Initial_pre_file"],"init_profile":inputs["StartProfilesData"][target.target]["StartProfilesFile"],"cwr":target.add["ComputeWithRadiation"],"td":target.add["Thermodiffusion"],"etf":target.add["ExpTempFile"],"copy":inputs["StartProfilesData"][target.target]["CopyTo"],"fuel":inputs["StartProfilesData"][target.target]["fuel"],"oxidizer":inputs["StartProfilesData"][target.target]["oxidizer"],"bathGas":inputs["StartProfilesData"][target.target]["bathGas"],"glob_rxn":inputs["StartProfilesData"][target.target]["globalReaction"],"to_Pa":target.pressure,"to_K":target.temperature,"units":yaml.dump(yaml.load(str(inputs["StartProfilesData"][target.target]["units"]),Loader=Loader),default_flow_style=True),"flow":target.add["flow_rate"],"id":target.uniqueID}

	infile = open("profile_generator.input",'w')
	infile.write(instring)
	infile.close()
	
	run_file = open("run_profile",'w')
	s = """#!/bin/bash
python3 %(profile_generator_src)s profile_generator.input &> profile.log
"""% {"profile_generator_src":inputs["Bin"]["bin"]+"/profileGenerator.py"}
#&> Flame.log ; gnome-terminal && tail -f flame.log
	run_file.write(s)
	run_file.close()
	subprocess.call(["chmod","+x",'run_profile'])

def create_inputs_to_start(inputs,phi,fuel,oxidizer,bathgas,globeRxn,P,T,controls):
	instring = ""
	if inputs["exp_type"]["exp"] == "Fls":
		instring = """
############
# Numerics #
############

#### Newton solver ####
TimeDepFlag = {time_flag}
DeltaTStart = 1.0e-8
DeltaTMax = 1.0e5
UseNumericalJac is TRUE
UseSecondOrdJac is TRUE
UseModifiedNewton = TRUE

DampFlag = TRUE
LambdaMin = 1.0e-2

MaxIter = 5000
TolRes = 1.0e-15
TolDy = 1e-4

#### grid ####

DeltaNewGrid = 25
OneSolutionOneGrid = TRUE
initialgridpoints = {init_grid}
maxgridpoints = {max_grid}
q = -0.25
R = 60

########################
# Sensitivity Analysis #
########################

#ReactionFluxAnal is TRUE

#######
# I/O #
#######

#WriteEverySolution = TRUE
#PrintMolarFractions is TRUE
#AdditionalOutput is TRUE

OutputPath is ./output
StartProfilesFile is {s_p_loc}

#############
# Chemistry #
#############

MechanismFile is {p_f}
globalReaction is {g_r};

{fuel_is}
oxidizer is {oxidizer}

#########
# Flame #
#########

Flame is {flame}
ExactBackward is TRUE

{phi}

pressure = {pressure}

ComputeWithRadiation is {radiation_tag}
Thermodiffusion is {thermodiffusion_tag}

#######################
# Boundary conditions #
#######################

#ConstMassFlux is TRUE
#MassFlux = 0.3

Unburnt Side {{
	dirichlet {{
		t = {temperature}
		{dirichlet}
	}}
}}

{bounds}

""".format(time_flag=inputs["TimeFlag"],init_grid=inputs["initialgridpoints"],max_grid=inputs["maxgridpoints"],s_p_loc=inputs["StartProfilesFile"],p_f=inputs["MechanismFile"],g_r=globeRxn,fuel_is=fuel["fuelIs"],oxidizer=oxidizer["oxidizerIs"],flame=inputs["Flame"],phi = phi,pressure=P,radiation_tag=inputs["ComputeWithRadiation"],thermodiffusion_tag=inputs["Thermodiffusion"],temperature=T,dirichlet=controls["dirichlet"],bounds=controls["conc_bounds"])

	elif inputs["exp_type"]["exp"] == "Flf":
	
		instring = """############
# Numerics #
############

#### Newton solver ####

UseNumericalJac is TRUE
UseSecondOrdJac is TRUE
UseModifiedNewton = TRUE

DampFlag = TRUE
LambdaMin = 1.0e-2

TimeDepFlag is {time_flag}
DeltaTStart = 1e-5

MaxIter = 5000
TolRes = 1.0e-15
TolDy = 1e-4

#### grid ####

DeltaNewGrid = 25
OneSolutionOneGrid = TRUE
initialgridpoints = {init_grid}
maxgridpoints = {max_grid}
q = -0.25
R = 60

########################
# Sensitivity Analysis #
########################

#ReactionFluxAnal is TRUE

#######
# I/O #
#######
#AdditionalOutput is TRUE
WriteEverySolution is TRUE
#PrintMolarFractions is TRUE

OutputPath is ./output
StartProfilesFile is {s_p_loc}

#############
# Chemistry #
#############

MechanismFile is {p_f}
globalReaction is {g_r};

{fuel_is}
oxidizer is {oxidizer}

#########
# Flame #
#########

Flame is {flame}
ExactBackward is TRUE

#phi = 1.9272

pressure = {pressure}

ComputeWithRadiation is {radiation_tag}
Thermodiffusion is {thermodiffusion_tag}
TransModel is MonoAtomic

ExpTempFile is {etp}

#######################
# Boundary conditions #
#######################

ConstLewisNumber is TRUE
ConstMassFlux is TRUE

MassFlowRate = {flow_rate}

Unburnt Side {{
	dirichlet {{
		t = {temperature}
		{dirichlet}
	}}
}}

{bounds}
""".format(time_flag=inputs["TimeFlag"],init_grid=inputs["initialgridpoints"],max_grid=inputs["maxgridpoints"],s_p_loc=inputs["StartProfilesFile"],p_f=inputs["MechanismFile"],g_r=globeRxn,fuel_is=fuel["fuelIs"],oxidizer=oxidizer["oxidizerIs"],flame=inputs["Flame"],pressure=P,radiation_tag=inputs["ComputeWithRadiation"],thermodiffusion_tag=inputs["Thermodiffusion"],etp=inputs["ExpTempFile"],flow_rate=inputs["flowRate"],temperature=T,dirichlet=controls["dirichlet"],bounds=controls["conc_bounds"])

	infile = open("FlameMaster.input",'w')
	infile.write(instring)
	infile.close()
	#print(target.uniqueID)
	run_file = open("run_generate",'w')
	s = """#!/bin/bash
{Bin}/FlameMan &> Flame.log
""".format(Bin=inputs["bin"])
#&> Flame.log ; gnome-terminal && tail -f flame.log
	run_file.write(s)
	run_file.close()
	subprocess.call(["chmod","+x",'run_generate'])


def create_input_file(case,opt_dict,target,mech_file=None):
	#print(target.uniqueID)
	thermo_file_location = opt_dict["Locations"]["thermo_file"]
	trans_file_location = opt_dict["Locations"]["trans_file"]
	startProfile_location = opt_dict["StartProfilesData"]
	file_specific_command = "-f chemkin"
	if mech_file == None:
		mech_file = "mechanism.yaml"
	
	global_reaction = opt_dict["Inputs"]["global_reaction"]
	instring = ''
	fuel = target.fuel
	P = target.pressure
	#print(target.target)
	#for i in g_reaction:
	#	global_reaction += i

	# Defining the default variables with solver specific keys:
	# In later version this moduel is to be shifted in combustion_target_class
	# creating an input dictionary, lateron the dictionary will be part of target object
	
	#print(target.target)
	#print(target.uniqueID)
	if target.input_file != None:
		instring = open(target.input_file,'r').read()
			
	elif "JSR" in target.target:
		#print(target.add)
		instring ="""import pandas as pd
import time
import cantera as ct

gas = ct.Solution("{mech}")
reactor_temperature = {temperature}  # Kelvin
reactor_pressure = {pressure}  # in atm. This equals 1.06 bars
inlet_concentrations = {species_conc}
gas.TPX = reactor_temperature, reactor_pressure, inlet_concentrations
species = "{species}"
# Reactor parameters
residence_time = {residenceTime}  # s
reactor_volume = 1  # m3
max_simulation_time = {maxSimTime}
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
tau_file.write("#T(K)	mole fraction"+"\\n"+"{{}}	{{}}".format(reactor_temperature,MF))
tau_file.close()""".format(mech = mech_file,temperature = target.temperature,pressure=target.pressure,species_conc = target.species_dict,species=str(target.add["species"]),residenceTime = target.add["residenceTime"],maxSimTime = target.add["maxSimulationTime"])
		extract = """
		"""
	elif "Flf" in target.target:
		if "cantera" in target.add["solver"]:
			instring = """#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import time
import pandas as pd

gas = ct.Solution('{mech}')

reactorTemperature = {temperature} #Kelvin
reactorPressure = {pressure} #in atm. This equals 1.06 bars
concentrations = {species_conc}
gas.TPX = reactorTemperature, reactorPressure, concentrations 

residenceTime = {residenceTime} #s

reactorVolume = 1 #m3

pressureValveCoefficient = 0.01
maxPressureRiseAllowed = 0.01

maxSimulationTime = {maxSimTime} # seconds

fuelAirMixtureTank = ct.Reservoir(gas)
exhaust = ct.Reservoir(gas)

stirredReactor = ct.IdealGasReactor(gas, energy='off', volume=reactorVolume)

massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
										   downstream=stirredReactor,
										   mdot=stirredReactor.mass/residenceTime)

pressureRegulator = ct.Valve(upstream=stirredReactor,
							 downstream=exhaust,
							 K=pressureValveCoefficient)

reactorNetwork = ct.ReactorNet([stirredReactor])
columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
columnNames = ['pressure'] + columnNames

# use the above list to create a DataFrame
timeHistory = pd.DataFrame(columns=columnNames)
tic = time.time()
t = 0
while t < maxSimulationTime:
	t = reactorNetwork.step()
	
state = np.hstack([stirredReactor.thermo.P, 
				   stirredReactor.mass, 
				   stirredReactor.volume, 
				   stirredReactor.T, 
				   stirredReactor.thermo.X])

toc = time.time()

concentrations = stirredReactor.thermo.X

# Store the result in the dataframe that indexes by temperature
timeHistory.loc[reactorTemperature] = state
Flf = str(timeHistory['{species}']).split("\\n")[0].split("	")[1]
flf_file = open("output/flf.out",'w')
flf_file.write("#T(K)	flf(conc.)"+"\\n"+"{{}}	{{}}".format(reactorTemperature,Flf))
flf_file.close()""".format(mech = mech_file,temperature = target.temperature,pressure=target.pressure,species_conc = target.species_dict,species=target.add["species"],residenceTime = target.add["residenceTime"],maxSimTime = target.add["maxSimulationTime"])
			extract = """
		"""
	elif "RCM" in target.target:
		if "cantera" in target.add["solver"] and "non_reactive" not in target.add["RCM_type"]: 
			instring="""import cantera as ct
import numpy as np
from pyked import ChemKED
from urllib.request import urlopen
import yaml,time
import pandas as pd
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

def ignitionDelay(df,pList, species,cond="max",specified_value ="None;None",exp_conc = "None"):
	if cond == "max":
		T = df["temperature"]
		time = df.index.to_numpy()
		conc = df[species].to_numpy()
		valid_indices = np.where((T>criteria))		
		dt = np.diff(time)
		dT = np.diff(T.to_numpy())
		slope = dT/dt
		index = int(slope.argmax())
		peaks,_ = find_peaks(conc[valid_indices])
		if len(peaks)>1:
			peak_ind = abs(np.asarray(peaks) - index).argmin()
		elif len(peaks) == 1:
			peak_ind = 0
		else:
			peaks,_ = find_peaks(conc)
			valid_indices = [[0]]
			if len(peaks) == 1:
				peak_ind = 0
			else:
				raise AssertionError("No ignition!!")
		tau = time[valid_indices[0][0]-1+peaks[peak_ind]]
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

T_initial = {temperature}
P_initial = {pressure}
X_initial = {species_conc}

global criteria
if T_initial> 2000:
	criteria = T_initial + 5
else:
	criteria = T_initial + 50


gas = ct.Solution('{mech}')
gas.TPX = T_initial, P_initial, X_initial

reac = ct.IdealGasReactor(gas)
env = ct.Reservoir(ct.Solution('{mech}'))

df = pd.read_csv("{volume_profile}")
exp_time = df["time(s)"]
exp_volume = df["volume(cm3)"]


stateVariableNames = [reac.component_name(item) for item in range(reac.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)
#exp_time = dp.volume_history.time.magnitude
#try:
#	exp_volume = dp.volume_history.volume.magnitude
#except:
#	exp_volume = dp.volume_history.quantity.magnitude
	
keywords = {{"vproTime": exp_time, "vproVol": exp_volume}}
ct.Wall(reac, env, velocity=VolumeProfile(keywords));

netw = ct.ReactorNet([reac])
netw.max_time_step = np.min(np.diff(exp_time))

t = 0
pressureList = []
counter = 1;
while(gas.T < 1800):
	t = netw.step()
	if (counter%1 == 0):	   	
		timeHistory.loc[t] = netw.get_state().astype('float64')
		pressureList.append(gas.P)
	counter+=1

tic = time.time() 
#tau = ignitionDelay(timeHistory,pressureList,"{delay_def}","{delay_cond}","{specific_cond}",{exp_conc})
tok = time.time()
#Actual defination for RCM using temperature profile
tempe = timeHistory["temperature"].to_numpy()
time_ = timeHistory.index.to_numpy()
index1 = tempe.argmax()
t1 = time_[index1]
tempe_split = np.array_split(tempe,3)
tempe_2 = np.concatenate((tempe_split[0],tempe_split[1]))
index2 = tempe_2.argmax()
t2 = time_[index2]
tau = t1-t2
print(f"t1:{{time_[index1]}},{{tempe[index1]}}\\nt2:{{time_[index2]}},{{tempe_2[index2]}}")
fig = plt.figure()
plt.plot(time_[index1],tempe[index1],"o")
plt.plot(time_[index2],tempe_2[index2],"o")
plt.plot(time_,tempe,"-")
plt.savefig("Ignition_def.pdf")
tau_file = open("output/RCM.out",'w')
tau_file.write("#T(K)	tau(us)"+"\\n"+"{{}}	{{}}".format(T_initial,tau*(10**6)))
tau_file.close()""".format(mech = mech_file,temperature = target.temperature_i,pressure=target.pressure_i,species_conc = target.species_dict,exp_conc = target.add["exp_conc"][float(target.temperature)],delay_def = target.add["ign_delay_def"],volume_profile=target.add["volume_profile"][float(target.temperature)],delay_cond = target.add["ign_cond"],specific_cond = target.add["specific_cond"],saveAll=target.add["saveAll"])
			extract = """
		"""
		elif "CHEMKIN_PRO" in target.add["solver"]:
			instring = """#!/usr/bin/env python
# encoding: utf-8

import os
import csv
import numpy as np
from chemkin import ChemkinJob, getIgnitionDelay, getIgnitionDelayOH
import pandas as pd
from scipy.interpolate import CubicSpline, Akima1DInterpolator
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

def find_first_slope_transition(time, temperature, temp_min, temp_max):
	# Create a DataFrame for easier manipulation
	df = pd.DataFrame({{"time": time, "temperature": temperature}})

	# Filter the data to the specified temperature range
	filtered_df = df[(df["temperature"] >= temp_min) & (df["temperature"] <= temp_max)].reset_index(drop=True)

	# Compute the first derivative (slope) of the temperature
	filtered_df["slope"] = np.gradient(filtered_df["temperature"], filtered_df["time"])

	# Find the first slope transition from positive to a smaller positive value
	for i in range(1, len(filtered_df) - 1):
		if filtered_df["slope"].iloc[i - 1] > 0 and filtered_df["slope"].iloc[i] < filtered_df["slope"].iloc[i - 1] and filtered_df["slope"].iloc[i] > 0:
			# This is the first transition from positive to smaller positive slope
			return {{
				"time": filtered_df["time"].iloc[i],
				"temperature": filtered_df["temperature"].iloc[i],
				"slope": filtered_df["slope"].iloc[i]
			}}

	# If no transition is found, return None
	return None


def find_local_peak_in_range(time, temperature, time_min, time_max):
	# Create a DataFrame for easier manipulation
	df = pd.DataFrame({{"time": time, "temperature": temperature}})

	# Filter the data to the specified time range
	filtered_df = df[(df["time"] >= time_min) & (df["time"] <= time_max)].reset_index(drop=True)

	# Compute the first derivative (slope) of the temperature
	filtered_df["slope"] = np.gradient(filtered_df["temperature"], filtered_df["time"])

	# Find the local peak by detecting where the slope changes from positive to negative
	for i in range(1, len(filtered_df) - 1):
		if filtered_df["slope"].iloc[i - 1] > 0 and filtered_df["slope"].iloc[i + 1] < 0:
			# This is a local peak
			return {{
				"time": filtered_df["time"].iloc[i],
				"temperature": filtered_df["temperature"].iloc[i],
				"slope": filtered_df["slope"].iloc[i]
			}}

	# If no peak is found, return None
	return None
def find_critical_point(time, temperature, time_min, time_max):
   
	# Create a DataFrame for easier manipulation
	df = pd.DataFrame({{"time": time, "temperature": temperature}})

	# Filter the data to the specified time range
	filtered_df = df[(df["time"] >= time_min) & (df["time"] <= time_max)].reset_index(drop=True)

	# Compute the first derivative (slope)
	filtered_df["slope"] = np.gradient(filtered_df["temperature"], filtered_df["time"])
	
	# Compute the second derivative (rate of change of slope)
	filtered_df["slope_change"] = np.gradient(filtered_df["slope"], filtered_df["time"])

	# 1. Find the first point where the slope is decreasing to zero or becomes negative
	first_slope_decrease = None
	for i in range(1, len(filtered_df)):
		if filtered_df["slope"].iloc[i] < 0 and filtered_df["slope"].iloc[i-1] > 0:
			first_slope_decrease = {{
				"type": "slope_decrease",
				"time": filtered_df["time"].iloc[i],
				"temperature": filtered_df["temperature"].iloc[i],
				"slope": filtered_df["slope"].iloc[i]
			}}
			break

	# 2. Find the first local peak (local maximum) in the temperature profile
	peak_idx = filtered_df["temperature"].idxmax()  # Index of max temperature
	peak_point = {{
		"type": "peak",
		"time": filtered_df["time"].iloc[peak_idx],
		"temperature": filtered_df["temperature"].iloc[peak_idx],
		"slope": filtered_df["slope"].iloc[peak_idx]
	}}

	# Return the first detected significant point (slope decrease or peak)
	if first_slope_decrease:
		return first_slope_decrease
	else:
		return peak_point


def find_slope_change_points(time, temperature, time_min, time_max):
	
	# Create a DataFrame for easier manipulation
	df = pd.DataFrame({{"time": time, "temperature": temperature}})

	# Filter the data within the specified time range
	filtered_df = df[(df["time"] >= time_min) & (df["time"] <= time_max)].reset_index(drop=True)

	# Calculate the first derivative (slope)
	filtered_df["slope"] = np.gradient(filtered_df["temperature"], filtered_df["time"])

	# Find the first slope decrease to near zero
	first_slope_decrease_point = None
	for i in range(1, len(filtered_df)):
		if filtered_df["slope"].iloc[i] < filtered_df["slope"].iloc[i - 1] and filtered_df["slope"].iloc[i] > 0:
			first_slope_decrease_point = {{
				"time": filtered_df["time"].iloc[i],
				"temperature": filtered_df["temperature"].iloc[i],
				"slope": filtered_df["slope"].iloc[i]
			}}
			break

	# Find the first major slope change
	second_derivative = np.gradient(filtered_df["slope"], filtered_df["time"])
	major_slope_change_idx = np.argmax(np.abs(second_derivative))  # Maximum absolute second derivative
	first_major_slope_change_point = {{
		"time": filtered_df["time"].iloc[major_slope_change_idx],
		"temperature": filtered_df["temperature"].iloc[major_slope_change_idx],
		"slope": filtered_df["slope"].iloc[major_slope_change_idx]
	}}

	return {{
		"first_slope_decrease_point": first_slope_decrease_point,
		"first_major_slope_change_point": first_major_slope_change_point
	}}

def find_nearest_temperature(time_list, temperature_list, target_time):
	
	# Convert to numpy arrays for easier computation
	time_array = np.array(time_list)
	temperature_array = np.array(temperature_list)
	
	# Find the index of the nearest time value
	idx = np.abs(time_array - target_time).argmin()
	
	# Get the nearest time and its corresponding temperature
	nearest_time = time_array[idx]
	nearest_temperature = temperature_array[idx]
	
	return {{"time": nearest_time, "temperature": nearest_temperature}}

def find_peak_point(time, temperature):
	# Create a DataFrame for easier manipulation
	df = pd.DataFrame({{"time": time, "temperature": temperature}})
	
	# Compute the derivative (slope) using finite differences
	df["slope"] = np.gradient(df["temperature"], df["time"])
	
	# Find the index of the maximum temperature
	max_temp_idx = df["temperature"].idxmax()
	
	# Get the time and temperature at the peak
	peak_time = df.loc[max_temp_idx, "time"]
	peak_temperature = df.loc[max_temp_idx, "temperature"]
	
	return {{"time": peak_time, "temperature": peak_temperature, "slope": df.loc[max_temp_idx, "slope"]}}


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
		T = df["temperature"]
		time = df.index.to_numpy()
		conc = df[species].to_numpy()
		valid_indices = np.where((T>criteria))		
		dt = np.diff(time)
		dT = np.diff(T.to_numpy())
		slope = dT/dt
		index = int(slope.argmax())
		peaks,_ = find_peaks(conc[valid_indices])
		if len(peaks)>1:
			peak_ind = abs(np.asarray(peaks) - index).argmin()
		elif len(peaks) == 1:
			peak_ind = 0
		else:
			peaks,_ = find_peaks(conc)
			valid_indices = [[0]]
			if len(peaks) == 1:
				peak_ind = 0
			else:
				raise AssertionError("No ignition!!")
		tau = time[valid_indices[0][0]-1+peaks[peak_ind]]
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


def convert_and_transpose_ckcsv(input_filepath, output_filepath):
	with open(input_filepath, 'r') as ckcsv_file:
		reader = csv.reader(ckcsv_file)

		# Skip the first row (metadata) and read the actual data
		rows = list(reader)
		data_rows = rows[1:-2]  # Skip metadata row
		
		# Extract variables, units, and time series data
		variables = [row[0] for row in data_rows]	 # First column as variable names
		units = [row[1] for row in data_rows]		 # Second column as units
		data_values = [row[2:] for row in data_rows]  # Time series data from the third column onward
		headers = [f"{{row[0]}} {{row[1]}}" for row in data_rows]
		#print(data_values[-1])	
		# Transpose the data values so each variable has its data in columns
		transposed_data = list(zip(*data_values))
		with open(output_filepath, 'w', newline='') as csv_file:
			writer = csv.writer(csv_file)

			# Write headers (variable names)
			writer.writerow(headers)

			# Write units row
			#writer.writerow(units)

			# Write transposed data rows
			writer.writerows(transposed_data)

	print(f"Conversion and transposition complete. Data saved to {{output_filepath}}")

#mechanism in yaml format: {mech}
target_Temperature = {target_Temperature}
# Parameters
reactorTemperature = {temperature_K}  # Replace with your desired temperature in Kelvin
reactorPressure = {pressure}  # Replace with your desired pressure in pa
spec = {species_conc}
conc_X = list(spec.items()) # Example species concentrations, adjust as needed
global criteria
if reactorTemperature>= 1800:
	criteria = reactorTemperature + 7
elif reactorTemperature>= 1400 and reactorTemperature< 1800:
	criteria = reactorTemperature + 50
else:
	criteria = reactorTemperature + 100
# Set up the Chemkin simulation
currentDir = os.path.dirname(__file__).strip("\\n")
chemFile = os.path.join(currentDir, 'mechanism.inp')  # Update with the mechanism file path
tempDir = os.path.join(currentDir, 'output')

estimatedIgnitionDelayTime = {estimateTIG}
variableVolumeProfileType = "{VP_FILE_TYPE}"
if variableVolumeProfileType == "ck_file":
	VP = "{volume_profile}"
elif variableVolumeProfileType == "Dict" or variableVolumeProfileType == "csv_file":
	df = pd.read_csv("{volume_profile}")
	exp_time = df["time(s)"]
	exp_volume = df["volume(cm3)"]
	keywords = {{"vproTime": exp_time, "vproVol": exp_volume}}
	VP = VolumeProfile(keywords)
else:
	VP = "{volume_profile}"

# Initialize Chemkin job
job = ChemkinJob(
	name=conc_X[0][0],
	chemFile=chemFile,
	tempDir=tempDir,
)

job.preprocess()

# Write Chemkin input for a single simulation
input_file = job.writeInputHomogeneousBatch(
	problemType='constrainVandSolveE',  # Problem type in Chemkin
	reactants=conc_X,
	temperature=reactorTemperature,
	pressure=reactorPressure,
	endTime= {estimateTIG},  # Simulation end time in seconds
	variableVolume = True,
	variableVolumeProfile = VP,
	variableVolumeProfileType = "{VP_FILE_TYPE}"
)

# Run the Chemkin simulation
job.run(input_file, model='CKReactorGenericClosed', pro=True)
job.postprocess(sens=False, rop=False, all=True, transpose=False)

#tau = ignitionDelay(timeHistory,pressureList,"{delay_def}","{delay_cond}","{specific_cond}",{exp_conc})

# Calculate ignition delay
try:
	tau = getIgnitionDelayOH(job.ckcsvFile)
except ValueError:
	print("Value Error")
	tau = None

saveAll = {saveAll}
if saveAll:
	convert_and_transpose_ckcsv(job.ckcsvFile, "time_history.csv")

timeHistory = pd.read_csv("time_history.csv")
pressureList = timeHistory["Pressure  (bar)"]
#tau = ignitionDelay(timeHistory,pressureList,"CHV","max","None",)

#Actual defination for RCM using temperature profile
tempe = timeHistory["Temperature  (K)"].to_numpy()
time_ = timeHistory.index.to_numpy()
smooth_temperature = gaussian_filter1d(tempe, sigma=5)  # Adjust sigma for smoothness

# Plot
plt.plot(time_, tempe, label="Original")
plt.plot(time_, smooth_temperature, label="Smoothed")
plt.legend()
plt.savefig("Gaussian_smoothness.pdf")
#tempe = smooth_temperature
slope_change_point_1 = find_first_slope_transition(time_,smooth_temperature,0.8*target_Temperature,target_Temperature+10)#, direction='forward')
t1 = slope_change_point_1['time']
T_1 = slope_change_point_1['temperature']
slope_change_point_2 = find_peak_point(time_,tempe)
t2 = slope_change_point_2["time"]
T_2 = slope_change_point_2['temperature']
max_T_1 = find_slope_change_points(time_,tempe,t1,0.9*t2) #Find the maxima in between t1 and t2
if max_T_1['first_slope_decrease_point']:
	t1_dash = max_T_1['first_slope_decrease_point']['time']
	T_1_dash = max_T_1['first_slope_decrease_point']['temperature']
else:   
	print(t1)
	new_T_1 = find_local_peak_in_range(time_,tempe,t1,1.06*t1)
	t1_dash = new_T_1['time']
	T_1_dash = new_T_1['temperature']
max_T_2 = find_slope_change_points(time_,smooth_temperature,t2,1.5*t2)

if max_T_2['first_slope_decrease_point']:
	t2_dash = max_T_2['first_slope_decrease_point']['time']
	T_2_dash = max_T_2['first_slope_decrease_point']['temperature']
else:   
	new_T_2 = find_critical_point(time_,smooth_temperature,t2,1.5*t2)
	
	t2_dash = new_T_2['time']
	T_2_dash = new_T_2['temperature']


max_T_3 = find_nearest_temperature(time_,tempe,(t2+t2_dash)/2)
t3 = max_T_3['time']
T_3 = max_T_3['temperature']
print(f"EOC:\\t{{t1_dash}},{{T_1_dash}}")
tau = t3-t1_dash
fig = plt.figure()
plt.plot(t1_dash,T_1_dash,"o")
plt.plot(t3,T_3,"o")
plt.plot(time_,tempe,"-")
plt.savefig("Ignition_def.pdf")
# Output the results
if tau:
	print(f"Ignition delay at {{reactorTemperature}} K and {{reactorPressure}} atm: {{tau*1e3}} ms")
	with open("output/RCM.out", 'w') as tau_file:
		tau_file.write(f"#T(K)	tau(us)\\n{{T_1_dash}}	{{tau * 1e6}}\\n")

import glob

def delete_files(directory):
	# Patterns for files to delete
	patterns = ["*.out", "*.asc", "*.zip"]
	exception_file = "tau.out"
	
	for pattern in patterns:
		# Get all files matching the pattern in the directory
		files = glob.glob(os.path.join(directory, pattern))
		
		for file in files:
			# Skip the exception file
			if os.path.basename(file).lower() == exception_file.lower():
				continue
			
			# Delete the file
			try:
				os.remove(file)
				print(f"Deleted: {{file}}")
			except Exception as e:
				print(f"Error deleting {{file}}: {{e}}")

# Specify the directory containing the files
directory_path = "output/"
delete_files(directory_path)
""".format(mech = mech_file,target_Temperature = target.temperature,temperature_K = target.temperature_i,pressure=target.pressure_i,species_conc = target.species_dict,exp_conc = target.add["exp_conc"][float(target.temperature)], volume_profile=target.add["volume_profile"][float(target.temperature)],VP_FILE_TYPE = target.add["volume_profile_type"],estimateTIG = float(target.add["estimateTIG"]),delay_def = target.add["ign_delay_def"],delay_cond = target.add["ign_cond"],specific_cond = target.add["specific_cond"],saveAll=target.add["saveAll"])

		elif "cantera" in target.add["solver"] and "non_reactive" in target.add["RCM_type"]: 
				instring="""#!/usr/bin/python
from __future__ import division
from __future__ import print_function
import pandas as pd
import numpy as np
import time
import cantera as ct
from scipy.interpolate import CubicSpline, Akima1DInterpolator
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.signal import find_peaks
class TemperatureFromPressure(object):
	def __init__(self, pressure, T_initial, chem_file='species.cti', cti_source=None,init_X=None):
		if cti_source is None:
			gas = ct.Solution(chem_file)
		else:
			gas = ct.Solution(source=cti_source)
		if init_X is None:
			gas.TP = T_initial, pressure[0]*1e5
		else:
			gas.TPX = T_initial, pressure[0]*1e5,init_X
		initial_entropy = gas.entropy_mass
		self.temperature = np.zeros((len(pressure)))
		for i, p in enumerate(pressure):
			gas.SP = initial_entropy, p*1e5
			self.temperature[i] = gas.T

class VolumeProfile(object):
	
	def __init__(self, time, volume):
		# The time and volume are stored as lists in the keywords
		# dictionary. The volume is normalized by the first volume
		# element so that a unit area can be used to calculate the
		# velocity.
		self.time = np.array(time)
		self.volume = np.array(volume)/volume[0]

		# The velocity is calculated by the forward difference.
		# numpy.diff returns an array one element smaller than the
		# input array, so we append a zero to match the length of the
		# self.time array.
		self.velocity = np.diff(self.volume)/np.diff(self.time)
		self.velocity = np.append(self.velocity, 0)

	def __call__(self, t):
	   

		if t < self.time[-1]:
			# prev_time_point is the previous value in the time array
			# after the current simulation time
			prev_time_point = self.time[self.time <= t][-1]
			# index is the index of the time array where
			# prev_time_point occurs
			index = np.where(self.time == prev_time_point)[0][0]
			return self.velocity[index]
		else:
			return 0

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
		T = df["temperature"]
		time = df.index.to_numpy()
		conc = df[species].to_numpy()
		valid_indices = np.where((T>criteria))		
		dt = np.diff(time)
		dT = np.diff(T.to_numpy())
		slope = dT/dt
		index = int(slope.argmax())
		peaks,_ = find_peaks(conc[valid_indices])
		if len(peaks)>1:
			peak_ind = abs(np.asarray(peaks) - index).argmin()
		elif len(peaks) == 1:
			peak_ind = 0
		else:
			peaks,_ = find_peaks(conc)
			valid_indices = [[0]]
			if len(peaks) == 1:
				peak_ind = 0
			else:
				raise AssertionError("No ignition!!")
		tau = time[valid_indices[0][0]-1+peaks[peak_ind]]
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
	
def doNonReactive(gas,inp_time,inp_vol):
	gas.set_multiplier(0)
	#r = ct.Reactor(contents=gas)
	r = ct.IdealGasReactor(gas)
	reactorNetwork = ct.ReactorNet([r])
	env = ct.Reservoir(ct.Solution('air.yaml'))
	ct.Wall(r, env, A=1.0, velocity=VolumeProfile(inp_time, inp_vol))
	
	stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
	timeHistory = pd.DataFrame(columns=stateVariableNames)

	#t0 = time.time()
	estimatedIgnitionDelayTime = 1
	t = 0
	pressureList = []
	volumeList = []
	counter = 1;
	while(t < estimatedIgnitionDelayTime):
		t = reactorNetwork.step()
		#print(t)
		if (counter%1 == 0):
			timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
			pressureList.append(gas.P)
			volumeList.append(r.volume)
		counter+=1
	
	return timeHistory,pressureList,volumeList

VP_file_type = "{VP_FILE_TYPE}"
if VP_file_type == "ck_file":
	file_ = "{volumeProfile}"
	df = pd.read_csv(file_, delim_whitespace=True, header=None)
	df.columns = ["Tag", "Time", "Volume"]
	inp_time = df["Time"]
	inp_vol = df["Volume"]

else:
	file_ = "{volumeProfile}"
	df = pd.read_csv(file_, sep=",")
	inp_time = df["time(s)"]
	inp_vol = df["volume(cm3)"]

TDC_time = inp_time[list(inp_vol).index(min(inp_vol))]
time = []
temperature = []
pressure = []
#input_volume = volume
simulated_volume = []
end_temp = 2500.
end_time = 0.2

gas = ct.Solution('{mech}')
gas_nr = ct.Solution("{mech}")

reactorTemperature = {temperature} #Kelvin
reactorPressure = {pressure}
global criteria
if reactorTemperature>= 1800:
	criteria = reactorTemperature + 7
elif reactorTemperature>= 1400 and reactorTemperature< 1800:
	criteria = reactorTemperature + 50
else:
	criteria = reactorTemperature + 100
gas.TPX = reactorTemperature, reactorPressure,{species_conc}
gas_nr.TPX = reactorTemperature, reactorPressure,{species_conc}

timeHistory_nonreactive,pressureList_nonreactive,vol_nonreactive = doNonReactive(gas_nr,inp_time,inp_vol)


fig = plt.figure()
plt.plot(timeHistory_nonreactive,pressureList_nonreactive,"r-",label="Non-reactive pressure profile")
plt.savefig("p_profile.png",bbox_inches="tight")""".format(mech = mech_file,volumeProfile = target.add["volume_profile"][float(target.temperature)],temperature = target.temperature,pressure=target.pressure,species_conc = target.species_dict,exp_conc = target.add["exp_conc"][float(target.temperature)],VP_FILE_TYPE = target.add["volume_profile_type"],delay_def = target.add["ign_delay_def"],delay_cond = target.add["ign_cond"],specific_cond = target.add["specific_cond"],saveAll=target.add["saveAll"])
		
		else:
			raise AssertionError("Invalid solver name for RCM simulations !!!")
	elif "Tig" in target.target:
		#print(target.ignition_type)
		#print(target.target)
		if target.ig_mode == "RCM" and target.simulation == "VTIM":
			if "cantera" in target.add["solver"]: 
				instring="""#!/usr/bin/python
from __future__ import division
from __future__ import print_function
import pandas as pd
import numpy as np
import time
import cantera as ct
from scipy.interpolate import CubicSpline, Akima1DInterpolator
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.signal import find_peaks
class TemperatureFromPressure(object):
	def __init__(self, pressure, T_initial, chem_file='species.cti', cti_source=None,init_X=None):
		if cti_source is None:
			gas = ct.Solution(chem_file)
		else:
			gas = ct.Solution(source=cti_source)
		if init_X is None:
			gas.TP = T_initial, pressure[0]*1e5
		else:
			gas.TPX = T_initial, pressure[0]*1e5,init_X
		initial_entropy = gas.entropy_mass
		self.temperature = np.zeros((len(pressure)))
		for i, p in enumerate(pressure):
			gas.SP = initial_entropy, p*1e5
			self.temperature[i] = gas.T

class VolumeProfile(object):
	
	def __init__(self, time, volume):
		# The time and volume are stored as lists in the keywords
		# dictionary. The volume is normalized by the first volume
		# element so that a unit area can be used to calculate the
		# velocity.
		self.time = np.array(time)
		self.volume = np.array(volume)/volume[0]

		# The velocity is calculated by the forward difference.
		# numpy.diff returns an array one element smaller than the
		# input array, so we append a zero to match the length of the
		# self.time array.
		self.velocity = np.diff(self.volume)/np.diff(self.time)
		self.velocity = np.append(self.velocity, 0)

	def __call__(self, t):
	   

		if t < self.time[-1]:
			# prev_time_point is the previous value in the time array
			# after the current simulation time
			prev_time_point = self.time[self.time <= t][-1]
			# index is the index of the time array where
			# prev_time_point occurs
			index = np.where(self.time == prev_time_point)[0][0]
			return self.velocity[index]
		else:
			return 0

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
		T = df["temperature"]
		time = df.index.to_numpy()
		conc = df[species].to_numpy()
		valid_indices = np.where((T>criteria))		
		dt = np.diff(time)
		dT = np.diff(T.to_numpy())
		slope = dT/dt
		index = int(slope.argmax())
		peaks,_ = find_peaks(conc[valid_indices])
		if len(peaks)>1:
			peak_ind = abs(np.asarray(peaks) - index).argmin()
		elif len(peaks) == 1:
			peak_ind = 0
		else:
			peaks,_ = find_peaks(conc)
			valid_indices = [[0]]
			if len(peaks) == 1:
				peak_ind = 0
			else:
				raise AssertionError("No ignition!!")
		tau = time[valid_indices[0][0]-1+peaks[peak_ind]]
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
	
def doNonReactive(gas,inp_time,inp_vol):
	gas.set_multiplier(0)
	#r = ct.Reactor(contents=gas)
	r = ct.IdealGasReactor(gas)
	reactorNetwork = ct.ReactorNet([r])
	env = ct.Reservoir(ct.Solution('air.yaml'))
	ct.Wall(r, env, A=1.0, velocity=VolumeProfile(inp_time, inp_vol))
	
	stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
	timeHistory = pd.DataFrame(columns=stateVariableNames)

	#t0 = time.time()
	estimatedIgnitionDelayTime = 1
	t = 0
	pressureList = []
	volumeList = []
	counter = 1;
	while(t < estimatedIgnitionDelayTime):
		t = reactorNetwork.step()
		#print(t)
		if (counter%1 == 0):
			timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
			pressureList.append(gas.P)
			volumeList.append(r.volume)
		counter+=1
	
	return timeHistory,pressureList,volumeList


volume = np.genfromtxt("{volumeProfile}", delimiter=',')
inp_time = volume[:, 0]
inp_vol = volume[:, 1]
TDC_time = inp_time[list(inp_vol).index(min(inp_vol))]
time = []
temperature = []
pressure = []
input_volume = volume
simulated_volume = []
end_temp = 2500.
end_time = 0.2

gas = ct.Solution('{mech}')
gas_nr = ct.Solution("{mech}")

reactorTemperature = {temperature} #Kelvin
reactorPressure = {pressure}
global criteria
if reactorTemperature>= 1800:
	criteria = reactorTemperature + 7
elif reactorTemperature>= 1400 and reactorTemperature< 1800:
	criteria = reactorTemperature + 50
else:
	criteria = reactorTemperature + 100
gas.TPX = reactorTemperature, reactorPressure,{species_conc}
gas_nr.TPX = reactorTemperature, reactorPressure,{species_conc}

#timeHistory_nonreactive,pressureList_nonreactive,vol_nonreactive = doNonReactive(gas_nr,inp_time,inp_vol)

r = ct.IdealGasReactor(gas)
reactorNetwork = ct.ReactorNet([r])
env = ct.Reservoir(ct.Solution('air.yaml'))
ct.Wall(r, env, A=1.0, velocity=VolumeProfile(inp_time, inp_vol))

stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)

#t0 = time.time()
estimatedIgnitionDelayTime = 1.0
t = 0
pressureList = []
volumeList = []
counter = 1;
while(t < estimatedIgnitionDelayTime):
	t = reactorNetwork.step()
	#print(t)
	if (counter%1 == 0):
		timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
		pressureList.append(gas.P)
		volumeList.append(r.volume)
	counter+=1

tau = ignitionDelay(timeHistory,pressureList,"{delay_def}","{delay_cond}","{specific_cond}",{exp_conc})
#TDC_point = vol_nonreactive.index(min(vol_nonreactive))
#Tc =TemperatureFromPressure(pressureList_nonreactive,reactorTemperature,chem_file="mechanism.yaml",init_X={species_conc}).temperature[TDC_point]
#Pc = pressureList_nonreactive[TDC_point]

tau_file = open("output/tau.out",'w')
tau_file.write("#T(K)	tau(us)"+"\\n"+"{{}}	{{}}".format(reactorTemperature,(tau-TDC_time)*(10**6)))
tau_file.close()
#T_file = open("output/TP.out",'w')
#T_file.write("#T(K)	P(bar)"+"\\n"+"{{}}	{{}}".format(Tc,Pc/1e5))
#T_file.close()
{saveAll}timeHistory.to_csv("time_history.csv")""".format(mech = mech_file,volumeProfile = target.add["volumeProfile"],temperature = target.temperature,pressure=target.pressure,species_conc = target.species_dict,exp_conc = target.add["exp_conc"][float(target.temperature)],delay_def = target.add["ign_delay_def"],delay_cond = target.add["ign_cond"],specific_cond = target.add["specific_cond"],saveAll=target.add["saveAll"])
		elif target.ignition_type == "reflected":
			#print(target.add["BoundaryLayer"])
			if "cantera" in target.add["solver"] and target.add["BoundaryLayer"]== True:
				instring ="""#!/usr/bin/python
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
dt = {dt}
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
		T = df["temperature"]
		time = df.index.to_numpy()
		conc = df[species].to_numpy()
		valid_indices = np.where((T>criteria))		
		dt = np.diff(time)
		dT = np.diff(T.to_numpy())
		slope = dT/dt
		index = int(slope.argmax())
		peaks,_ = find_peaks(conc[valid_indices])
		if len(peaks)>1:
			peak_ind = abs(np.asarray(peaks) - index).argmin()
		elif len(peaks) == 1:
			peak_ind = 0
		else:
			peaks,_ = find_peaks(conc)
			valid_indices = [[0]]
			if len(peaks) == 1:
				peak_ind = 0
			else:
				raise AssertionError("No ignition!!")
		tau = time[valid_indices[0][0]-1+peaks[peak_ind]]
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
	time =np.arange(t1,t2,{dt})
	return time

	
def getPressureProfile(gas,time,dpdt):
	p = gas.P
	p_new = p
	PressureProfile = []
	dt = {dt}
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
gas = ct.Solution('{mech}')

reactorTemperature = {temperature} #Kelvin
reactorPressure = {pressure}

global criteria
if reactorTemperature>= 1800:
	criteria = reactorTemperature + 7
elif reactorTemperature>= 1400 and reactorTemperature< 1800:
	criteria = reactorTemperature + 50
else:
	criteria = reactorTemperature + 100

gas.TPX = reactorTemperature, reactorPressure,{species_conc}
r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
env = ct.Reservoir(ct.Solution('{mech}'))

dpdt = {dpdt}#%/ms
estimatedIgnitionDelayTime = {estimateTIG}
flag = False
time = getTimeProfile(0,estimatedIgnitionDelayTime)
pressure_profile = getPressureProfile(gas,time,dpdt)
volume_profile = getVolumeProfile_From_PressureProfile(gas,pressure_profile)
string = "time(s),volume(cm3)\\n"
for i, t in enumerate(time):
	string+=f"{{t}},{{volume_profile[i]}}\\n"
g = open(f"VTIM_P_{{int(reactorPressure/100000)}}_T_{{int(reactorTemperature)}}.csv","w").write(string)

keywords = {{"vproTime": time, "vproVol": volume_profile}}
ct.Wall(r, env, velocity=VolumeProfile(keywords));

reactorNetwork = ct.ReactorNet([r])

reactorNetwork.max_time_step = 0.0001

stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)


####============================#####
### Solver specific information #####
###	Do not change			#####
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
#	t = reactorNetwork.step()
#	#print(t,gas.T)
#	time_profile.append(t)
#	gas_T.append(gas.T)
#	time_diff = np.diff(np.asarray(time_profile))
#	conc = np.diff(np.asarray(gas_T))
#	if len(conc) == 1:
#		slope = conc/time_diff
#		index = slope[0]
#	else:
#		slope = conc/time_diff
#		index = int(slope.argmax())
#	slope_arg_max.append(index)
#	count+=1
#	if max(slope_arg_max) == slope_arg_max[-1] and count==1:
#		continue
#	elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature <= 200 and t<estimatedIgnitionDelayTime:
#		continue
#	elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature <= 200 and t>estimatedIgnitionDelayTime:
#		break
#	elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature >= 200 and t>estimatedIgnitionDelayTime:
#		break
#	elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature >= 100 and t<estimatedIgnitionDelayTime and flag == False:
#		estimatedIgnitionDelayTime = t + 1e-4
#		flag = True
#	if (counter%1 == 0):
#		#if reactorNetwork.get_state().astype('float64')[index_a] <0:
#		 #   raise AssertionError("Numerical error!!")			
#		timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
#		pressureList.append(gas.P)
#	counter+=1
tau = ignitionDelay(timeHistory,pressureList,"{delay_def}","{delay_cond}","{specific_cond}",{exp_conc})
# Toc
tau_file = open("output/tau.out",'w')
tau_file.write("#T(K)	tau(us)"+"\\n"+"{{}}	{{}}".format(reactorTemperature,tau*(10**6)))
tau_file.close()
species = "{delay_def}"
time = timeHistory.index.to_numpy()
if species != "p":
	conc = (timeHistory[species]).to_numpy()
else:
	conc = pressureList
	
fig = plt.figure()
plt.plot(time,conc,"b-",label="{delay_def} profile")
plt.legend()
plt.savefig("profile.pdf")
saveAll = {saveAll}
if saveAll == True:
	timeHistory.to_csv("time_history.csv")
""".format(mech = mech_file,temperature = target.temperature,pressure=target.pressure,species_conc = target.species_dict,exp_conc = target.add["exp_conc"][float(target.temperature)],dpdt = target.add["dpdt"],dt = float(target.add["dt"]),estimateTIG = float(target.add["estimateTIG"]), delay_def = target.add["ign_delay_def"],delay_cond = target.add["ign_cond"],specific_cond = target.add["specific_cond"],saveAll=target.add["saveAll"])
				
			elif "cantera" in target.add["solver"] and target.add["BoundaryLayer"] != True:
				instring = """#!/usr/bin/python
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

gas = ct.Solution('{mech}')

reactorTemperature = {temperature} #Kelvin
reactorPressure = {pressure}

global criteria
if reactorTemperature>= 1800:
	criteria = reactorTemperature + 7
elif reactorTemperature>= 1400 and reactorTemperature< 1800:
	criteria = reactorTemperature + 50
else:
	criteria = reactorTemperature + 100

gas.TPX = reactorTemperature, reactorPressure,{species_conc}
r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
reactorNetwork = ct.ReactorNet([r])
reactorNetwork.max_time_step = 0.0001
stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)
#index_a = stateVariableNames.index("{delay_def}")

def ignitionDelay(df,pList, species,cond="max",specified_value ="None;None",exp_conc = "None"):
	if cond == "max":
		T = df["temperature"]
		time = df.index.to_numpy()
		conc = df[species].to_numpy()
		valid_indices = np.where((T>criteria))		
		dt = np.diff(time)
		dT = np.diff(T.to_numpy())
		slope = dT/dt
		index = int(slope.argmax())
		peaks,_ = find_peaks(conc[valid_indices])
		if len(peaks)>1:
			peak_ind = abs(np.asarray(peaks) - index).argmin()
		elif len(peaks) == 1:
			peak_ind = 0
		else:
			peaks,_ = find_peaks(conc)
			valid_indices = [[0]]
			if len(peaks) == 1:
				peak_ind = 0
			else:
				raise AssertionError("No ignition!!")
		tau = time[valid_indices[0][0]-1+peaks[peak_ind]]
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
estimatedIgnitionDelayTime = {estimateTIG}
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
#	t = reactorNetwork.step()
#	#print(t,gas.T)
#	time_profile.append(t)
#	gas_T.append(gas.T)
#	time_diff = np.diff(np.asarray(time_profile))
#	conc = np.diff(np.asarray(gas_T))
#	if len(conc) == 1:
#		slope = conc/time_diff
#		index = slope[0]
#	else:
#		slope = conc/time_diff
#		index = int(slope.argmax())
#	slope_arg_max.append(index)
#	count+=1
#	if max(slope_arg_max) == slope_arg_max[-1] and count==1:
#		continue
#	elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature <= 200 and t<estimatedIgnitionDelayTime:
#		continue
#	elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature <= 200 and t>estimatedIgnitionDelayTime:
#		break
#	elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature >= 200 and t>estimatedIgnitionDelayTime:
#		break
#	elif max(slope_arg_max) == slope_arg_max[-1] and max(slope_arg_max) == slope_arg_max[-2] and gas.T - reactorTemperature >= 100 and t<estimatedIgnitionDelayTime and flag == False:
#		estimatedIgnitionDelayTime = t + 1e-4
#		flag = True
#	if (counter%1 == 0):
#		#if reactorNetwork.get_state().astype('float64')[index_a] <0:
#		 #   raise AssertionError("Numerical error!!")			
#		timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
#		pressureList.append(gas.P)
#	counter+=1

tau = ignitionDelay(timeHistory,pressureList,"{delay_def}","{delay_cond}","{specific_cond}",{exp_conc})
# Toc
t1 = time.time()
tau_file = open("output/tau.out",'w')
tau_file.write("#T(K)	tau(us)"+"\\n"+"{{}}	{{}}".format(reactorTemperature,tau*(10**6)))
tau_file.close()
saveAll = {saveAll}
if saveAll == True:
	timeHistory.to_csv("time_history.csv")
""".format(mech = mech_file,temperature = target.temperature,pressure=target.pressure,species_conc = target.species_dict,exp_conc = target.add["exp_conc"][float(target.temperature)],delay_def = target.add["ign_delay_def"],delay_cond = target.add["ign_cond"],estimateTIG = target.add["estimateTIG"],specific_cond = target.add["specific_cond"],saveAll=target.add["saveAll"])
			
			elif "CHEMKIN_PRO" in target.add["solver"] and target.add["BoundaryLayer"] != True:
				
				
				instring = """#!/usr/bin/env python
# encoding: utf-8

import os
import csv
import numpy as np
from chemkin import ChemkinJob, getIgnitionDelay, getIgnitionDelayOH
import pandas as pd
from scipy.interpolate import CubicSpline, Akima1DInterpolator
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

def ignitionDelay(df,pList, species,cond="max",specified_value ="None;None",exp_conc = "None"):
	if cond == "max":
		T = df["temperature"]
		time = df.index.to_numpy()
		conc = df[species].to_numpy()
		valid_indices = np.where((T>criteria))		
		dt = np.diff(time)
		dT = np.diff(T.to_numpy())
		slope = dT/dt
		index = int(slope.argmax())
		peaks,_ = find_peaks(conc[valid_indices])
		if len(peaks)>1:
			peak_ind = abs(np.asarray(peaks) - index).argmin()
		elif len(peaks) == 1:
			peak_ind = 0
		else:
			peaks,_ = find_peaks(conc)
			valid_indices = [[0]]
			if len(peaks) == 1:
				peak_ind = 0
			else:
				raise AssertionError("No ignition!!")
		tau = time[valid_indices[0][0]-1+peaks[peak_ind]]
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


def convert_and_transpose_ckcsv(input_filepath, output_filepath):
	with open(input_filepath, 'r') as ckcsv_file:
		reader = csv.reader(ckcsv_file)

		# Skip the first row (metadata) and read the actual data
		rows = list(reader)
		data_rows = rows[1:-2]  # Skip metadata row
		
		# Extract variables, units, and time series data
		variables = [row[0] for row in data_rows]	 # First column as variable names
		units = [row[1] for row in data_rows]		 # Second column as units
		data_values = [row[2:] for row in data_rows]  # Time series data from the third column onward
		headers = [f"{{row[0]}} {{row[1]}}" for row in data_rows]
		#print(data_values[-1])	
		# Transpose the data values so each variable has its data in columns
		transposed_data = list(zip(*data_values))
		with open(output_filepath, 'w', newline='') as csv_file:
			writer = csv.writer(csv_file)

			# Write headers (variable names)
			writer.writerow(headers)

			# Write units row
			#writer.writerow(units)

			# Write transposed data rows
			writer.writerows(transposed_data)

	print(f"Conversion and transposition complete. Data saved to {{output_filepath}}")

#mechanism in yaml format: {mech}

# Parameters
reactorTemperature = {temperature_K}  # Replace with your desired temperature in Kelvin
reactorPressure = {pressure}  # Replace with your desired pressure in atm
spec = {species_conc}
conc_X = list(spec.items()) # Example species concentrations, adjust as needed

# Set up the Chemkin simulation
currentDir = os.path.dirname(__file__).strip("\\n")
chemFile = os.path.join(currentDir, 'mechanism.inp')  # Update with the mechanism file path
tempDir = os.path.join(currentDir, 'output')
global criteria
if reactorTemperature>= 1800:
	criteria = reactorTemperature + 7
elif reactorTemperature>= 1400 and reactorTemperature< 1800:
	criteria = reactorTemperature + 50
else:
	criteria = reactorTemperature + 100

# Initialize Chemkin job
job = ChemkinJob(
	name=conc_X[0][0],
	chemFile=chemFile,
	tempDir=tempDir,
)

job.preprocess()

# Write Chemkin input for a single simulation
input_file = job.writeInputHomogeneousBatch(
	problemType='constrainVandSolveE',  # Problem type in Chemkin
	reactants=conc_X,
	temperature=reactorTemperature,
	pressure=reactorPressure,
	endTime= {End_Time}  # Simulation end time in seconds
)

# Run the Chemkin simulation
job.run(input_file, model='CKReactorGenericClosed', pro=True)
job.postprocess(sens=False, rop=False, all=True, transpose=False)

# Calculate ignition delay
#try:
#	tau = getIgnitionDelayOH(job.ckcsvFile)
#except ValueError:
#	print("Value Error")
#	tau = None

saveAll = True 
if saveAll:
	convert_and_transpose_ckcsv(job.ckcsvFile, "time_history.csv")

timeHistory = pd.read_csv("time_history.csv",index_col=0)
pressureList = timeHistory["Pressure  (bar)"]
species = f"Mole_fraction_{delay_def}  ()"
tau = ignitionDelay(timeHistory,pressureList,species,"{delay_cond}","{specific_cond}",{exp_conc})

# Output the results
if tau:
	print(f"Ignition delay at {{reactorTemperature}} K and {{reactorPressure}} atm: {{tau*1e3}} ms")
	with open("output/tau.out", 'w') as tau_file:
		tau_file.write(f"#T(K)	tau(us)\\n{{reactorTemperature}}	{{tau * 1e6}}\\n")

import glob

def delete_files(directory):
	# Patterns for files to delete
	patterns = ["*.out", "*.asc", "*.zip"]
	exception_file = "tau.out"
	
	for pattern in patterns:
		# Get all files matching the pattern in the directory
		files = glob.glob(os.path.join(directory, pattern))
		
		for file in files:
			# Skip the exception file
			if os.path.basename(file).lower() == exception_file.lower():
				continue
			
			# Delete the file
			try:
				os.remove(file)
				print(f"Deleted: {{file}}")
			except Exception as e:
				print(f"Error deleting {{file}}: {{e}}")

# Specify the directory containing the files
directory_path = "output/"
delete_files(directory_path)
""".format(mech = mech_file,temperature_K = target.temperature,pressure=target.pressure,species_conc = target.species_dict,exp_conc = target.add["exp_conc"][float(target.temperature)],End_Time = target.add["estimateTIG"], delay_def = target.add["ign_delay_def"],delay_cond = target.add["ign_cond"],specific_cond = target.add["specific_cond"],saveAll=target.add["saveAll"])
				

			elif "CHEMKIN_PRO" in target.add["solver"] and target.add["BoundaryLayer"] == True:
				instring = """#!/usr/bin/env python
# encoding: utf-8

import os
import csv
import numpy as np
from chemkin import ChemkinJob, getIgnitionDelay, getIgnitionDelayOH
import pandas as pd
from scipy.interpolate import CubicSpline, Akima1DInterpolator
import cantera as ct
from scipy.signal import find_peaks

global dt
dt = {dt}
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
		T = df["temperature"]
		time = df.index.to_numpy()
		conc = df[species].to_numpy()
		valid_indices = np.where((T>criteria))		
		dt = np.diff(time)
		dT = np.diff(T.to_numpy())
		slope = dT/dt
		index = int(slope.argmax())
		peaks,_ = find_peaks(conc[valid_indices])
		if len(peaks)>1:
			peak_ind = abs(np.asarray(peaks) - index).argmin()
		elif len(peaks) == 1:
			peak_ind = 0
		else:
			peaks,_ = find_peaks(conc)
			valid_indices = [[0]]
			if len(peaks) == 1:
				peak_ind = 0
			else:
				raise AssertionError("No ignition!!")
		tau = time[valid_indices[0][0]-1+peaks[peak_ind]]
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
	time =np.arange(t1,t2,{dt})
	return time

	
def getPressureProfile(gas,time,dpdt):
	p = gas.P
	p_new = p
	PressureProfile = []
	dt = {dt}
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


def convert_and_transpose_ckcsv(input_filepath, output_filepath):
	with open(input_filepath, 'r') as ckcsv_file:
		reader = csv.reader(ckcsv_file)

		# Skip the first row (metadata) and read the actual data
		rows = list(reader)
		data_rows = rows[1:-2]  # Skip metadata row
		
		# Extract variables, units, and time series data
		variables = [row[0] for row in data_rows]	 # First column as variable names
		units = [row[1] for row in data_rows]		 # Second column as units
		data_values = [row[2:] for row in data_rows]  # Time series data from the third column onward
		headers = [f"{{row[0]}} {{row[1]}}" for row in data_rows]
		#print(data_values[-1])	
		# Transpose the data values so each variable has its data in columns
		transposed_data = list(zip(*data_values))
		with open(output_filepath, 'w', newline='') as csv_file:
			writer = csv.writer(csv_file)

			# Write headers (variable names)
			writer.writerow(headers)

			# Write units row
			#writer.writerow(units)

			# Write transposed data rows
			writer.writerows(transposed_data)

	print(f"Conversion and transposition complete. Data saved to {{output_filepath}}")

#mechanism in yaml format: {mech}
gas = ct.Solution('{mech}')

# Parameters
reactorTemperature = {temperature_K}  # Replace with your desired temperature in Kelvin
reactorPressure = {pressure}  # Replace with your desired pressure in atm
spec = {species_conc}
conc_X = list(spec.items()) # Example species concentrations, adjust as needed
gas.TPX = reactorTemperature, reactorPressure,spec
# Set up the Chemkin simulation
currentDir = os.path.dirname(__file__).strip("\\n")
chemFile = os.path.join(currentDir, 'mechanism.inp')  # Update with the mechanism file path
tempDir = os.path.join(currentDir, 'output')

global criteria
if reactorTemperature>= 1800:
	criteria = reactorTemperature + 7
elif reactorTemperature>= 1400 and reactorTemperature< 1800:
	criteria = reactorTemperature + 50
else:
	criteria = reactorTemperature + 100

dpdt = {dpdt}#%/ms
estimatedIgnitionDelayTime = {estimateTIG}
time = getTimeProfile(0,estimatedIgnitionDelayTime)
pressure_profile = getPressureProfile(gas,time,dpdt)
volume_profile = getVolumeProfile_From_PressureProfile(gas,pressure_profile)

keywords = {{"vproTime": time, "vproVol": volume_profile}}
VP = VolumeProfile(keywords)

# Initialize Chemkin job
job = ChemkinJob(
	name=conc_X[0][0],
	chemFile=chemFile,
	tempDir=tempDir,
)

job.preprocess()

# Write Chemkin input for a single simulation
input_file = job.writeInputHomogeneousBatch(
	problemType='constrainVandSolveE',  # Problem type in Chemkin
	reactants=conc_X,
	temperature=reactorTemperature,
	pressure=reactorPressure,
	endTime= estimatedIgnitionDelayTime,  # Simulation end time in seconds
	variableVolume = True,
	variableVolumeProfile = VP,
	variableVolumeProfileType = "Dict"
)

# Run the Chemkin simulation
job.run(input_file, model='CKReactorGenericClosed', pro=True)
job.postprocess(sens=False, rop=False, all=True, transpose=False)


# Calculate ignition delay
#try:
#	tau = getIgnitionDelayOH(job.ckcsvFile)
#except ValueError:
#	print("Value Error")
#	tau = None


saveAll = True 
if saveAll:
	convert_and_transpose_ckcsv(job.ckcsvFile, "time_history.csv")

timeHistory = pd.read_csv("time_history.csv",index_col=0)
pressureList = timeHistory["Pressure  (bar)"]
species = f"Mole_fraction_{delay_def}  ()"
tau = ignitionDelay(timeHistory,pressureList,species,"{delay_cond}","{specific_cond}",{exp_conc})

# Output the results
if tau:
	print(f"Ignition delay at {{reactorTemperature}} K and {{reactorPressure}} atm: {{tau*1e3}} ms")
	with open("output/tau.out", 'w') as tau_file:
		tau_file.write(f"#T(K)	tau(us)\\n{{reactorTemperature}}	{{tau * 1e6}}\\n")


import glob

def delete_files(directory):
	# Patterns for files to delete
	patterns = ["*.out", "*.asc", "*.zip"]
	exception_file = "tau.out"
	
	for pattern in patterns:
		# Get all files matching the pattern in the directory
		files = glob.glob(os.path.join(directory, pattern))
		
		for file in files:
			# Skip the exception file
			if os.path.basename(file).lower() == exception_file.lower():
				continue
			
			# Delete the file
			try:
				os.remove(file)
				print(f"Deleted: {{file}}")
			except Exception as e:
				print(f"Error deleting {{file}}: {{e}}")

# Specify the directory containing the files
directory_path = "output/"
delete_files(directory_path)
""".format(mech = mech_file,temperature_K = target.temperature,pressure=target.pressure,species_conc = target.species_dict,exp_conc = target.add["exp_conc"][float(target.temperature)], dpdt = target.add["dpdt"],dt = float(target.add["dt"]),estimateTIG = float(target.add["estimateTIG"]),delay_def = target.add["ign_delay_def"],delay_cond = target.add["ign_cond"],specific_cond = target.add["specific_cond"],saveAll=target.add["saveAll"])

			elif "FlameMaster" in target.add["solver"] and target.add["BoundaryLayer"] == True:
				instring = """############
# Numerics #
############

RelTol = 1.0e-9
AbsTol = 1.0e-12

TStart = 0.0
TEnd = 0.002

########################
# Sensitivity Analysis #
########################

#SensAnalReac is TRUE
#SensAnalSpec is TRUE

#FirstSensRate = 5

#SensMax is TRUE
#SensFinal is TRUE

#SensObjAll is TRUE
#SensObj is OH
#SensObj is H

#SensAnalFac = 2.0


#######
# I/O #
#######

#AdditionalOutput is TRUE
WriteEverySolution is TRUE
PrintMolarFractions is TRUE

OutputPath is ./output
NOutputs = 2000

#############
# Chemistry #
#############

MechanismFile is mechanism.pre
{}
globalReaction is {g_r};

{fuel_is}
oxidizer is O2

#########
# Flame #
#########

Flame is {simulation}
ExactBackward is TRUE

#phi = {phi}

Pressure = {pressure}

PressureChange = {dpdt}

#######################
# Boundary conditions #
#######################

#ContInc = -25
#ContType is Temperature
#ContBound = 800

InitialCond {{
	t = {temperature}
	{init_condition}
	
}}
	""".format(pressure = target.pressure, init_condition = target.initialCond , temperature = target.temperature, phi = target.phi, simulation = target.simulation, fuel_is = target.fuel_is, g_r = global_reaction, dpdt = target.add["dpdt"][target.temperature]/100)	
				extract = """#!/usr/bin/python
import numpy as np
import os
import pandas as pd
import scipy as sp
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import CubicSpline, Akima1DInterpolator
import matplotlib.pyplot as plt
def find_nearest(array, value):
	array = np.asarray(np.abs(np.log(array)))
	value = np.abs(np.log(value))
	idx = (np.abs(array - value)).argmin()
	return idx,array[idx]
def find_nearest_value(xi,x,y):
	index = min(range(len(x)),key=lambda i: abs(x[i] - xi))
	return y[index]	
os.chdir("output")
list_files = os.listdir()

ignitionDelayDefination = "{ign_def}"
ign_cond = "{cond}"
specific = "{specific}"
atomic_weight = {molecular_weight}

def ignitionDelay(df,species,cond="max",specified_value ="None;None"):
	for i in df:
	#if "p" in ignitionDelayDefination:
		a = i.strip()
		if a == "X-"+species:
			Y = list(df[i])
			break
		if a == "t[ms]":
			X = list(df[i])
	time = np.asarray(X)
	profile = np.asarray(Y)
	t = time
	P = profile
	if cond == "max":
		tau = time[list(profile).index(max(profile))]
		
	elif cond == "onset":
		time = np.diff(time)
		conc = np.diff(profile)
		slope = conc/time
		index = int(np.diff(slope).argmax())
		tau = t[index]
	elif cond == "dt-max":
		time = np.diff(time)
		conc = np.diff(profile)
		slope = conc/time
		index = int(slope.argmax())
		tau = t[index]
	elif cond == "specific":
		if specified_value.split(";")[0] == None:
			raise Assertionerror("Input required for specified_value in ignition delay")
		else:
			target = float(specified_value.split(";")[0])
			unit = specified_value.split(";")[1]
			if unit == "molecule":
				avogadro = 6.02214E+23
				target = target/avogadro
				molecular_wt = atomic_weight
				target = molecular_wt*target
				#index = np.where(np.abs(profile-np.ones(len(profile))*target)<0.15*target)[0]
				#tau = time[index[0]]
				index,value = find_nearest(profile,target)
				tau = time[index]	
			else:
				#print(conc)
				f = CubicSpline(time,profile, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))
				
				#print(time_new)
				
				conc_new = f(time_new)
				
				index,value = find_nearest(conc_new,target)
				#print(conc_new[root_index])
				tau = time[index]
			
	return tau

for i in list_files:
	if "p" in ignitionDelayDefination:
		if i.startswith("Y"):
			outfile = i
			break
		else:
			continue
	else:
		if i.startswith("X"):
			outfile = i
			break
		else:
			continue
data = []
read = open(outfile,"r").readlines()
for line in read[1:]:
	data.append(line.strip("\\n"))
file_dump = open("modified.csv","w+")
file_dump.write(("\\n").join(data))
file_dump.close()
file_dump = open("../time_history.csv").write(("\\n").join(data))

df = pd.read_csv("modified.csv",sep="\\t")
for i in df:
	#print(i)
	if "p" in ignitionDelayDefination:
		if "P" in i:
			profile = df[i]
			continue
	if "X-"+ignitionDelayDefination+" " in i:
		profile = df[i]
		
	if "t" in i:
		time = df[i]
	
tau = ignitionDelay(df,ignitionDelayDefination,ign_cond,specific)
file_d2 = open("tau.out","w+")
file_d2.write("tau\\t{{}}\\tms".format(tau))
file_d2.close()
""".format(molecular_weight=target.add["mol_wt"],ign_def = target.add["ign_delay_def"],cond = target.add["ign_cond"],specific=target.add["specific_cond"])
			else:
				instring = """############
# Numerics #
############

RelTol = 1.0e-9
AbsTol = 1.0e-12

TStart = 0.0
TEnd = 0.002

########################
# Sensitivity Analysis #
########################

#SensAnalReac is TRUE
#SensAnalSpec is TRUE

#FirstSensRate = 5

#SensMax is TRUE
#SensFinal is TRUE

#SensObjAll is TRUE
#SensObj is OH
#SensObj is H

#SensAnalFac = 2.0


#######
# I/O #
#######

#AdditionalOutput is TRUE
WriteEverySolution is TRUE
PrintMolarFractions is TRUE

OutputPath is ./output
NOutputs = 2000

#############
# Chemistry #
#############

MechanismFile is mechanism.pre
globalReaction is {g_r};

{fuel_is}
oxidizer is O2

#########
# Flame #
#########

Flame is {simulation}
ExactBackward is TRUE

#phi = {phi}

Pressure = {pressure}


#######################
# Boundary conditions #
#######################

#ContInc = -25
#ContType is Temperature
#ContBound = 800

InitialCond {{
	t = {temperature}
	{init_condition}
	
}}
	""".format(pressure = target.pressure, init_condition = target.initialCond , temperature = target.temperature, phi = target.phi, simulation = target.simulation, fuel_is = target.fuel_is, g_r = global_reaction)
			extract = """#!/usr/bin/python
import numpy as np
import os
import pandas as pd
import scipy as sp
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import CubicSpline, Akima1DInterpolator
import matplotlib.pyplot as plt
def find_nearest(array, value):
	array = np.asarray(np.abs(np.log(array)))
	value = np.abs(np.log(value))
	idx = (np.abs(array - value)).argmin()
	return idx,array[idx]
def find_nearest_value(xi,x,y):
	index = min(range(len(x)),key=lambda i: abs(x[i] - xi))
	return y[index]   
os.chdir("output")
list_files = os.listdir()

ignitionDelayDefination = "{ign_def}"
ign_cond = "{cond}"
specific = "{specific}"
atomic_weight = {molecular_weight}

from scipy.signal import find_peaks
def ignitionDelay(df,species,cond="max",specified_value ="None;None"):
	for i in df:
		a = i.strip()
		if a == "X-"+species:
			Y = list(df[i])
		elif "P" in i:
			Y = list(df[i])
		elif "T[K]" in i:
			T = list(df[i])
		else:
			#continue
			print(f"Data for {{species}} not found")
		if a == "t[ms]":
			X = list(df[i])
	time = np.asarray(X)
	conc = np.asarray(Y)
	T = np.asarray(T)
	profile = np.asarray(Y)
	P = profile
	if cond == "max":
		#T = df["T[K]"]
		time = df.index.to_numpy()
		#conc = df[species].to_numpy()
		dt = np.diff(time)
		dT = np.diff(T)
		slope = dT/dt
		index = int(slope.argmax())
		peaks,_ = find_peaks(conc)
		if len(peaks)>1:
			peak_ind = abs(np.asarray(peaks) - index).argmin()
		elif len(peaks) == 1:
			peak_ind = 0
		else:
			raise AssertionError("No ignition!!")
		tau = time[peaks[peak_ind]]
		#tau = time[list(profile).index(max(profile))]
		
	elif cond == "onset" and species != "p":
		dt = np.diff(time)
		dconc = np.diff(conc)
		#print(dconc)
		slope = dconc/dt
		index = int(np.diff(slope).argmax())
		tau = time[index]
	elif cond == "onset" and species == "p":
		dt = np.diff(time)
		dP = np.diff(profile)
		slope = dP/dt
		index = int(np.diff(slope).argmax())
		tau = time[index]
	elif cond == "dt-max" and species != "p":
		dt = np.diff(time)
		dconc = np.diff(conc)
		slope = dconc/dt
		index = int(slope.argmax())
		tau = time[index]
	elif cond == "dt-max" and species == "p":
		dt = np.diff(time)
		dP = np.diff(profile)
		slope = dP/dt
		index = int(slope.argmax())
		tau = time[index]
	elif cond == "specific":
		if specified_value.split(";")[0] == None:
			raise Assertionerror("Input required for specified_value in ignition delay")
		else:
			target = float(specified_value.split(";")[0])
			unit = specified_value.split(";")[1]
			if unit == "molecule":
				avogadro = 6.02214E+23
				target = target/avogadro
				molecular_wt = atomic_weight
				target = molecular_wt*target
				#index = np.where(np.abs(profile-np.ones(len(profile))*target)<0.15*target)[0]
				#tau = time[index[0]]
				index,value = find_nearest(profile,target)
				tau = time[index]	
			elif unit == "molecule/cm3":
				unit_conversion = 10**(6)
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
				if exp_conc != "":
					target = target/exp_conc				
				f = Akima1DInterpolator(time,conc)
				time_new = np.arange(min(time),max(time),10**(-8))
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)			
	return tau
	
for i in list_files:
	if "p" in ignitionDelayDefination:
		if i.startswith("Y"):
			outfile = i
			break
		else:
			continue
	else:
		if i.startswith("X"):
			outfile = i
			break
		else:
			continue
data = []
read = open(outfile,"r").readlines()
for line in read[1:]:
	data.append(line.strip("\\n"))
file_dump = open("modified.csv","w+")
file_dump.write(("\\n").join(data))
file_dump.close()
file_dump = open("../time_history.csv",'w').write(("\\n").join(data))
df = pd.read_csv("modified.csv",sep="\\t")
for i in df:
	#print(i)
	if "p" in ignitionDelayDefination:
		if "P" in i:
			profile = df[i]
			continue
	if "X-"+ignitionDelayDefination+" " in i:
		profile = df[i]
		
	if "t" in i:
		time = df[i]
	
tau = ignitionDelay(df,ignitionDelayDefination,ign_cond,specific)
file_d2 = open("tau.out","w+")
file_d2.write("tau\\t{{}}\\tms".format(tau))
file_d2.close()
""".format(molecular_weight=target.add["mol_wt"],ign_def = target.add["ign_delay_def"],cond = target.add["ign_cond"],specific=target.add["specific_cond"])
		else:
			raise AssertionError("Invalid ignition mode specified for ignition delay simulations")
	elif "Fls" in target.target:
		if "cantera" in target.add["solver"] and "phi" not in target.add["type"]:
			instring = """#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import cantera as ct
import numpy as np
import pandas as pd
To 	= {temperature}
Po  = {pressure}
gas = ct.Solution('{mech}')
gas.TPX = To, Po,{species_conc}
#gas.set_equivalence_ratio(1.0, "CH4", {{"O2": 1.0, "N2": 3.76}})
width = {width}
flame = ct.FreeFlame(gas, width=width)
flame.set_refine_criteria(ratio={ratio}, slope={slope}, curve={curve}, prune = 0.003 )#3,0.015,0.015
loglevel = {loglevel}
refine_grid = True
flame.soret_enabled = True
flame.energy_enabled = True
#flame.max_grid_points = 800
flame.transport_model = "Multi"
flame.solve(loglevel=loglevel,refine_grid = True, auto={auto})#1,True
Su0 = flame.velocity[0] #m/s
Su_file = open("output/Su.out",'w')
Su_file.write("#T(K)	Su(cm/s)"+"\\n"+"{{}}	{{}}".format(To,Su0*(100)))
Su_file.close()
""".format(mech = mech_file,temperature = target.temperature,pressure=target.pressure,species_conc = target.species_dict,width = target.add["width"],ratio = target.add["ratio"],slope = target.add["slope"],curve = target.add["curve"],loglevel = target.add["loglevel"],auto = target.add["auto"])
				
		elif "cantera" in target.add["solver"] and "phi" in target.add["type"]:
			instring = """#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import cantera as ct
import numpy as np
import pandas as pd
To 	= {temperature}
Po  = {pressure}
gas = ct.Solution('{mech}')
#gas.TPX = To, Po,{species_conc}

gas.set_equivalence_ratio({phi}, {fuel}, {{"O2": 1.0, "N2": 3.76}})
width = {width}
flame = ct.FreeFlame(gas, width=width)
flame.set_refine_criteria(ratio={ratio}, slope={slope}, curve={curve}, prune = 0.003 )#3,0.015,0.015
loglevel = {loglevel}
refine_grid = True
flame.soret_enabled = True
flame.energy_enabled = True
#flame.max_grid_points = 800
flame.transport_model = "Multi"
flame.solve(loglevel=loglevel,refine_grid = True, auto={auto})#1,True
Su0 = flame.velocity[0] #m/s
Su_file = open("output/Su.out",'w')
Su_file.write("#T(K)	Su(cm/s)"+"\\n"+"{{}}	{{}}".format(To,Su0*(100)))
Su_file.close()
""".format(mech = mech_file,temperature = target.temperature,phi = target.phi,fuel= target.add["fuel"],pressure=target.pressure,species_conc = target.species_dict,width = target.add["width"],ratio = target.add["ratio"],slope = target.add["slope"],curve = target.add["curve"],loglevel = target.add["loglevel"],auto = target.add["auto"])
		
		else:
			if target.uniqueID not in os.listdir(startProfile_location[target.target]["CopyTo"]):
				opt_dict = copy.deepcopy(opt_dict)
				opt_dict["StartProfilesData"][target.target]["fuel"]["key"] = str(opt_dict["StartProfilesData"][target.target]["fuel"]["key"])+"-"+str(target.fuel_type)
				
				opt_dict["StartProfilesData"][target.target]["fuel"]["To"] = target.fuel_id
				opt_dict["StartProfilesData"][target.target]["fuel"]["ToConc"] = target.fuel_x
				opt_dict["StartProfilesData"][target.target]["oxidizer"]["To"] = target.oxidizer
				opt_dict["StartProfilesData"][target.target]["oxidizer"]["ToConc"] = target.oxidizer_x 
				opt_dict["StartProfilesData"][target.target]["bathGas"]["key"] =str(opt_dict["StartProfilesData"][target.target]["bathGas"]["key"])+"-"+str(target.bath_gas_id) 
				opt_dict["StartProfilesData"][target.target]["bathGas"]["To"] = target.bath_gas
				opt_dict["StartProfilesData"][target.target]["bathGas"]["ToConc"] = target.bath_gas_x
				opt_dict["StartProfilesData"][target.target]["globalReaction"]["ToRxn"] = opt_dict["Inputs"]["global_reaction"]
				opt_dict["StartProfilesData"][target.target]["fuel"] = yaml.dump(yaml.load(str(opt_dict["StartProfilesData"][target.target]["fuel"]),Loader=Loader),default_flow_style=True)
				opt_dict["StartProfilesData"][target.target]["oxidizer"] = yaml.dump(yaml.load(str(opt_dict["StartProfilesData"][target.target]["oxidizer"]),Loader=Loader),default_flow_style=True)
				opt_dict["StartProfilesData"][target.target]["bathGas"] = yaml.dump(yaml.load(str(opt_dict["StartProfilesData"][target.target]["bathGas"]),Loader=Loader),default_flow_style=True)
				opt_dict["StartProfilesData"][target.target]["globalReaction"] = yaml.dump(yaml.load(str(opt_dict["StartProfilesData"][target.target]["globalReaction"]),Loader=Loader),default_flow_style=True)
				target.add["Fsc"] = file_specific_command
				create_start_profile_input(opt_dict,target)
				os.mkdir("output")
				subprocess.call("./run_profile")
				target.add["StartProfilesFile"] = startProfile_location[target.target]["CopyTo"]+"/"+str(target.uniqueID)
				
				#raise AssertionError("NOOOO!!")
			else:
				target.add["StartProfilesFile"] = startProfile_location[target.target]["CopyTo"]+"/"+str(target.uniqueID)
			instring = """############
# Numerics #
############

#### Newton solver ####
DeltaTStart = 1.0e-8
DeltaTMax = 1.0e5
UseNumericalJac is TRUE
UseSecondOrdJac is TRUE
UseModifiedNewton = TRUE

DampFlag = TRUE
LambdaMin = 1.0e-2

MaxIter = 5000
TolRes = 1.0e-15
TolDy = 1e-4

#### grid ####

DeltaNewGrid = 25
OneSolutionOneGrid = TRUE
initialgridpoints = 89
maxgridpoints = 700
q = -0.25
R = 60

########################
# Sensitivity Analysis #
########################

#ReactionFluxAnal is TRUE

#######
# I/O #
#######

WriteEverySolution is TRUE
PrintMolarFractions is TRUE
AdditionalOutput is TRUE

OutputPath is ./output
StartProfilesFile is {s_p_loc}

#############
# Chemistry #
#############

MechanismFile is mechanism.pre
globalReaction is {g_r};

{fuel_is}
oxidizer is {oxidizer}

#########
# Flame #
#########

Flame is {simulation}
ExactBackward is TRUE

#{phi}

pressure = {pressure}

ComputeWithRadiation is {radiation_tag}
Thermodiffusion is {thermodiffusion_tag}

#######################
# Boundary conditions #
#######################

#ConstMassFlux is TRUE
#MassFlux = 0.3

Unburnt Side {{
	dirichlet {{
		t = {temperature}
		{dirichlet}
	}}
}}

""".format(pressure = target.pressure, dirichlet = target.initialCond , temperature = target.temperature, phi = target.phi, simulation = target.simulation, fuel_is = target.fuel_is,oxidizer=target.oxidizer, g_r = global_reaction,s_p_loc=target.add["StartProfilesFile"],radiation_tag=target.add["ComputeWithRadiation"],thermodiffusion_tag=target.add["Thermodiffusion"])
#"time_flag":inputs["TimeFlag"],"init_grid":inputs["initialgridpoints"],"max_grid":inputs["maxgridpoints"],"s_p_loc":inputs["StartProfilesFile"],"p_f":inputs["MechanismFile"],"g_r":globeRxn,"fuel_is":fuel["fuelIs"],"oxidizer":oxidizer["oxidizerIs"],"flame":inputs["Flame"],"phi":phi,"pressure":P,"radiation_tag":inputs["ComputeWithRadiation"],"thermodiffusion_tag":inputs["Thermodiffusion"],"temperature":T,"dirichlet":controls["dirichlet"],"conc_bounds":controls["conc_bounds"])
			
			extract = """
		"""
	
	elif "Flf_a" in target.target:
		if "cantera" in target.add["solver"]:
				
			instring = """#!/usr/bin/python
import cantera as ct
import numpy as np
from pathlib import Path
import pandas as pd
import os


def printSolution(df,species,criteria):
	for i in df:
		if species == i:
			if "max" in criteria:
				x = "{{}}({{}})".format(criteria,species)
				string = "{{}}\t{{}}".format(x,max(df[i]))
	return string

p = {pressure}
tburner = {temperature}
#To available from extrapolating the temperature profile
#One must provide the temperature profile, and using the extrapolation should 
# write a code to
# find the temperature at x = 0
# using extrapolation techniques
mdot = {mdot} #mass flux to the burner
reactants = {species_conc}  # premixed gas composition
width = {width_bsf} # cm
width = width/100 #m
loglevel = {loglevel}  # amount of diagnostic output (0 to 5)
refine_grid = True  # 'True' to enable refinement
gas = ct.Solution('{mech}')
gas.TPX = tburner, p, reactants
os.chdir('output')
f = ct.BurnerFlame(gas=gas, width=width)
f.burner.mdot = mdot
data_file = "{data_file}"
zloc, tvalues = np.genfromtxt(str(data_file), delimiter=',', comments='#').T
zloc /= max(zloc)
f.set_refine_criteria(ratio={ratio}, slope={slope_bsf}, curve={curve})
# set the temperature profile to the values read in
f.flame.set_fixed_temp_profile(zloc, tvalues)

# show the initial estimate for the solution
f.show_solution()

# don't solve the energy equation
f.energy_enabled = False

# first solve the flame with mixture-averaged transport properties
f.transport_model = '{transport_model}'
f.set_refine_criteria(ratio={ratio}, slope={slope_bsf}, curve={curve})

#f.solve(loglevel, refine_grid)
f.solve({solve_bsf})
try:
	# save to HDF container file if h5py is installed
	f.write_hdf('burner_flame.h5', group='{group}', mode='w',
				description='{description}')
except ImportError:
	f.save('burner_flame.xml', '{group}', '{description}')
f.write_csv('burner_flame.csv', quiet=False)
f.show_stats()
df = pd.read_csv("burner_flame.csv")
species = 'X_{target}'
criteria = "{criteria}"
string = printSolution(df,species,criteria)
file_dump = open("result.dout","w+")
file_dump.write(string)
file_dump.close()""".format(mech = mech_file,temperature = target.burner_temp,pressure=target.pressure,species_conc = target.species_dict,width_bsf = target.add["flf_grid"],ratio = target.add["ratio"],slope_bsf = target.add["slope_bsf"],curve = target.add["curve"],loglevel = target.add["loglevel"],solve_bsf = target.add["solve_bsf"],transport_model = target.add["transport_model"],group = target.add["group"],description = target.add["description"],data_file=target.add["ExpTempFile"],mdot=target.add["flow_rate"],target = target.add["flf_target"],criteria = target.add["flf_cond"])
					
		else:
			if target.uniqueID not in os.listdir(startProfile_location[target.target]["CopyTo"]):
				print("No startProfile found for target {}".format(target.uniqueID))
				#print(opt_dict["StartProfilesData"][target.target]["CopyTo"])
				print("---------------")
				print("Generating new startProfile")
				opt_dict["StartProfilesData"][target.target]["fuel"]["key"] = str(opt_dict["StartProfilesData"][target.target]["fuel"]["key"])+"-"+str(target.fuel_type)
				opt_dict["StartProfilesData"][target.target]["fuel"]["To"] = target.fuel_id
				opt_dict["StartProfilesData"][target.target]["fuel"]["ToConc"] = target.fuel_x
				opt_dict["StartProfilesData"][target.target]["oxidizer"]["To"] = target.oxidizer
				opt_dict["StartProfilesData"][target.target]["oxidizer"]["ToConc"] = target.oxidizer_x 
				opt_dict["StartProfilesData"][target.target]["bathGas"]["key"] =str(opt_dict["StartProfilesData"][target.target]["bathGas"]["key"])+"-"+str(target.bath_gas_id) 
				opt_dict["StartProfilesData"][target.target]["bathGas"]["To"] = target.bath_gas
				opt_dict["StartProfilesData"][target.target]["bathGas"]["ToConc"] = target.bath_gas_x
				opt_dict["StartProfilesData"][target.target]["globalReaction"]["ToRxn"] = opt_dict["Inputs"]["global_reaction"]
				opt_dict["StartProfilesData"][target.target]["fuel"] = yaml.dump(yaml.load(str(opt_dict["StartProfilesData"][target.target]["fuel"]),Loader=Loader),default_flow_style=True)
				opt_dict["StartProfilesData"][target.target]["oxidizer"] = yaml.dump(yaml.load(str(opt_dict["StartProfilesData"][target.target]["oxidizer"]),Loader=Loader),default_flow_style=True)
				opt_dict["StartProfilesData"][target.target]["bathGas"] = yaml.dump(yaml.load(str(opt_dict["StartProfilesData"][target.target]["bathGas"]),Loader=Loader),default_flow_style=True)
				opt_dict["StartProfilesData"][target.target]["globalReaction"] = yaml.dump(yaml.load(str(opt_dict["StartProfilesData"][target.target]["globalReaction"]),Loader=Loader),default_flow_style=True)
				target.add["Fsc"] = file_specific_command
				create_start_profile_input(opt_dict,target)
				os.mkdir("output")
				subprocess.call("./run_profile")
				target.add["StartProfilesFile"] = startProfile_location[target.target]["CopyTo"]+"/"+str(target.uniqueID)
			else:
				target.add["StartProfilesFile"] = startProfile_location[target.target]["CopyTo"]+"/"+str(target.uniqueID)
		
			instring = """############
# Numerics #
############

#### Newton solver ####

UseNumericalJac is TRUE
UseSecondOrdJac is TRUE
UseModifiedNewton = TRUE

DampFlag = TRUE
LambdaMin = 1.0e-2
DeltaTStart = 1e-5

MaxIter = 5000
TolRes = 1.0e-15
TolDy = 1e-4

#### grid ####

DeltaNewGrid = 25
OneSolutionOneGrid = TRUE
initialgridpoints = 89
maxgridpoints = 300
q = -0.25
R = 60

########################
# Sensitivity Analysis #
########################

#ReactionFluxAnal is TRUE

#######
# I/O #
#######
AdditionalOutput is TRUE
WriteEverySolution is TRUE
PrintMolarFractions is TRUE

OutputPath is ./output
StartProfilesFile is {s_p_loc}

#############
# Chemistry #
#############

MechanismFile is mechanism.pre
globalReaction is {g_r};

{fuel_is}
oxidizer is {oxidizer}

#########
# Flame #
#########

Flame is {simulation}
ExactBackward is TRUE

#phi = 1.9272

pressure = {pressure}

ComputeWithRadiation is {radiation_tag}
Thermodiffusion is {thermodiffusion_tag}
TransModel is MonoAtomic

ExpTempFile is {etp}

#######################
# Boundary conditions #
#######################

ConstLewisNumber is TRUE
ConstMassFlux is TRUE

MassFlowRate = {flow_rate}

Unburnt Side {{
	dirichlet {{
		t = {temperature}
		{dirichlet}
	}}
}}

""".format(pressure = target.pressure,oxidizer=target.oxidizer, dirichlet = target.initialCond , temperature = target.temperature, phi = target.phi, simulation = target.simulation, fuel_is = target.fuel_is, g_r = global_reaction,s_p_loc=target.add["StartProfilesFile"],thermodiffusion_tag=target.add["Thermodiffusion"],radiation_tag=target.add["ComputeWithRadiation"],flow_rate=target.add["flow_rate"],etp=target.add["ExpTempFile"])
			extract = """
		"""
	
	elif "Flw" in target.target:
	#No startProfile needed
		if target.reactor == "JSR" or target.reactor =="PSR" or target.reactor =="FlowReactor":
			if "cantera" in target.add["solver"]:
				
				instring = """#!/usr/bin/python
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json,os
########################################################################
## Input Parameters
########################################################################

T_0 = {temperature}  # inlet temperature [K]
pressure = {pressure} # constant pressure [Pa]
composition_0 = {species_conc}
length = {flowReactorLength}  # *approximate* PFR length [m]
#u_0 = {flow_velocity}  # inflow velocity [m/s]
area = {crossSectionalArea}  # cross-sectional area [m**2]
reactorVolume = {reactorVolume}
## input file containing the reaction mechanism
reaction_mechanism = {mech}
residenceTime = {residence_time}
species = "{species_in_investigation}"
## Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
## of 'n_steps' stirred reactors.
n_steps = {time_step}
######################################################################
#env = ct.Reservoir(ct.Solution('air.yaml'))
#####################################################################
# Method 1: Lagrangian Particle Simulation
# With heat loss to the surroundings through
# walls
#####################################################################
# A Lagrangian particle is considered which travels through the PFR. Its
# state change is computed by upwind time stepping. The PFR result is produced
# by transforming the temporal resolution into spatial locations.
# The spatial discretization is therefore not provided a priori but is instead
# a result of the transformation.

# import the gas model and set the initial conditions
gas1 = ct.Solution(reaction_mechanism)
gas1.TPX = T_0, pressure, composition_0
#mass_flow_rate1 = u_0 * gas1.density * area

# create a new reactor
r1 = ct.IdealGasConstPressureReactor(gas1)
#gas1, energy = "on",volume =reactorVolume
#mdot = {flow_rate}
mdot = r1.mass/residenceTime
# If iso-thermal is false
#w = ct.Wall(r1, env, A=1.0, U=0.05)

# create a reactor network for performing time integration
sim1 = ct.ReactorNet([r1])

# approximate a time step to achieve a similar resolution as in the next method
t_total = {total_time}
dt = t_total / n_steps
# define time, space, and other information vectors
t1 = (np.arange(n_steps) + 1) * dt
z1 = np.zeros_like(t1)
u1 = np.zeros_like(t1)
states1 = ct.SolutionArray(r1.thermo)
for n1, t_i in enumerate(t1):
	# perform time integration
	sim1.advance(t_i)
	# compute velocity and transform into space
	u1[n1] = mdot/ area / r1.thermo.density
	z1[n1] = z1[n1 - 1] + u1[n1] * dt
	states1.append(r1.thermo.state)
os.chdir("output")
#####################################################################
#plt.figure()
#plt.plot(t1, states1.X[:, gas1.species_index(species)], label='Lagrangian Particle')
#plt.xlabel('$t$ [s]')
#plt.ylabel('$X_{{}}$ [-]'.format(species))
#plt.legend(loc=0)
#plt.savefig('pfr_X{{}}_t.png'.format(species))


df_t = pd.DataFrame(t1)
df_t.to_csv("states_time.csv")
df_X = pd.DataFrame(states1.X[:, gas1.species_index(species)],columns = [species])
df_X.to_csv("states_X.csv")
df_T = pd.DataFrame(states1.T,columns=["T"])
df_T.to_csv("states_T.csv")
#Creating method to extract the slope of species decomposition
#create a def for interpolation:
time_scale = []
for i in t1:
	if i<=residenceTime:
		time_scale.append(i)
	elif i == residenceTime:
		time_scale.append(i)
	else:
		continue
time_range = np.asarray(time_scale)
 
def interpolation2D(d, x):
	output = d[0][1] + (x - d[0][0]) * ((d[1][1] - d[0][1])/(d[1][0] - d[0][0]))
	return output

def getNearest(value,array):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return idx
#create a def to populate the array:
def populate(target,array,x):
	new = []
	for i,ele in enumerate(array):
		if i<=x:
			new.append(ele)
			continue
		elif i==x+1:
			new.append(target)
			new.append(ele)
			continue
		else:
			new.append(ele)
			continue
	new = np.asarray(new)
	return new
#create a def to extract the details from the time array
def getTimeStamp(target,array):
	if target in array:
		for i in array:
			if i == target:
				array = list(array)
				x = array.index(i)
				y = len(array[x+1:])
				array = np.asarray(array)
	else:	   
		idx = getNearest(target,array)
		if array[idx] > target:
		#between idx-1 and idx
			x = idx-1
			y = len(array[x:])
			new = populate(target,array,x)
			array = new	
		#can use chebychev interpolation technique
		#start with linear interpolation
		else: 
		#array[idx] < target:
		#between idx and idx+1
			x = idx
			y = len(array[x:])
			new = populate(target,array,x)
			array = new
	return array,x,y

def timeScaleShift(target,array,x,y):
	if target in array:
		for i in array:
			if i == target:
				array = list(array)
				mid = array.index(i)
				#print(mid)
				middle = [array[mid]]
				top = array[mid-x:mid]
				bottom = array[mid+1:mid+y+1]
				scale_shift = np.asarray(top+middle+bottom)
	else:		 
		idx = getNearest(target,array)
		if array[idx] > target:
		#between idx-1 and idx
			x_ = idx-1
			y_ = idx
			new = populate(target,array,x_-1)
			array = new
			scale_shift = timeScaleShift(target,array,x,y)
#		 
#can use chebychev interpolation technique
#		 #start with linear interpolation
		else: 
#		 #array[idx] < target:
#		 #between idx and idx+1
			x_ = idx
			y_ = idx+1
			new = populate(target,array,x_-1)
			#print(new)
			array = new
			scale_shift = timeScaleShift(target,array,x,y)
			#			 array = new
	return scale_shift

def rateDefination(g1,g2):
	t1= g1[0]
	x1 =g1[1]
	t2 = g2[0]
	x2 = g2[1]
	rate = np.abs((x2-x1)/((t1-t2)*1000))
	dt = abs(t2-t1)*1000
	#rate calculated in ppm/ms
	#time calculated in ms
	return rate,dt

def getDataPoints(arr_t,arr_x,X):
	if X in arr_x:
		for i,x in enumerate(arr_x):
			if x == X:
				x_p = x
				t_p = arr_t[i]
	else:
		idx = getNearest(X,arr_x)
		if idx == len(arr_t):
			idx = idx-1
		if arr_x[idx] > X:
			x1_ = idx-1
			x2_ = idx
			t1 = arr_t[idx-1]
			t2 = arr_t[idx]
			t_p = interpolation2D(np.array([[arr_x[x1_],t1],[arr_x[x2_],t2]]),X)
			x_p = X  
#		 
#can use chebychev interpolation technique
#		 #start with linear interpolation
		else: 
#		 #array[idx] < target:
#		 #between idx and idx+1
			x1_ = idx
			x2_ = idx+1
			t1 = arr_t[idx]
			t2 = arr_t[idx+1]
			t_p = interpolation2D(np.array([[arr_x[x1_],t1],[arr_x[x2_],t2]]),X)
			x_p = X  
	return (t_p,x_p)   
	
def getRate(xo,t_half,X,t,add):
	if "slope" in add["method"]:	
		if "percentage" in add["unit"]:
			fact = 100
		else:
			fact = 1
		i_y1 = float(add["range_"][0])/fact
		i_y2 = float(add["range_"][1])/fact
		i_anchor = float(add["anchor"])/fact
		i_array_x = X
		i_array_t = t
		ni_array_t,len_x,len_y = getTimeStamp(t_half,i_array_t)
		ni_array_x = timeScaleShift(float(xo)*i_anchor,i_array_x,len_x,len_y)
		p1 = getDataPoints(i_array_t,ni_array_x,float(xo)*i_y1)
		p2 = getDataPoints(i_array_t,ni_array_x,float(xo)*i_y2)
		rate,dt = rateDefination(p1,p2)
	return ni_array_t,ni_array_x,rate,dt

def load(species,string):
	string = str(string).strip("{{}}")
	st = string.split(",")
	for i in st:
		k = i.split(":")
		if str(species) == k[0].strip('""').strip("''"):
			return_val = float(k[1])
	return return_val
	
xo = load(species,composition_0)
X = df_X.to_numpy().flatten()
t = time_range

t_half = {anchor_time}
add = dict(method = "{method}",range_ = {limits},anchor = {anchor},unit = "{unit}")
ni_array_t,ni_array_x,rate,dt = getRate(xo,t_half,X,t,add)

rate_ = open("rate.csv","w")
stringR = "rate	{{}}	ppm/ms".format(float(rate*1000000))
rate_.write(stringR)
rate_.close()
time_ = open("time.csv","w")
stringT = "time	{{}}	ms".format(dt)
time_.write(stringT)
time_.close()
#plt.figure()
#plt.plot(t,ni_array_x, label='timeScaleShift')
#plt.xlabel('$t$ [s]')
#plt.ylabel('$X_{{}}$ [-]'.format(species))
#plt.xlim(0,residenceTime)
#plt.legend(loc=0)
#plt.savefig('compare_exp.png')
os.chdir("..")""".format(mech = mech_file,temperature = target.temperature ,pressure=target.pressure,species_conc = target.species_dict ,flowReactorLength = target.add["flw_length"] ,flow_velocity = target.add["flow_velocity"] ,crossSectionalArea = target.add["crossSectionalArea"] ,reactorVolume = target.add["reactorVolume"] ,residence_time = target.add["residenceTime"] ,total_time = target.add["total_time"],time_step = target.add["time_step"],anchor_time = target.add["anchor_time"] ,flow_rate = target.flow_rate ,species_in_investigation = target.add["flw_species"] ,method = target.add["flw_method"] ,limits = target.add["flw_limits"] ,anchor = target.add["anchor"] ,unit = target.add["limit_units"])
				
			else:
				instring = """############
# Numerics #
############

RelTol = 1.0e-10
AbsTol = 1.0e-12

TStart = 0.0
TEnd = {tend}
TRes = {tres}
########################
# Sensitivity Analysis #
########################

#SensAnalReac is TRUE
#SensAnalSpec is TRUE

#FirstSensRate = 5

#SensMax is TRUE
#SensFinal is TRUE

#SensObjAll is TRUE
#SensObj is OH
#SensObj is H
SensObjALL is TRUE
#SensAnalFac = 2.0


#######
# I/O #
#######

#AdditionalOutput is TRUE
WriteEverySolution is TRUE
PrintMolarFractions is TRUE

OutputPath is ./output
NOutputs = 1000

#############
# Chemistry #
#############

MechanismFile is mechanism.pre
globalReaction is {g_r};

{fuel_is}
oxidizer is O2

#########
# Flame #
#########

Flame is {simulation} 

{isotherm}

{heat}

AmbientTemp is {T_amb}


ExactBackward is TRUE

#{phi}

Pressure = {pressure}


#######################
# Boundary conditions #
#######################

#ContInc = -25
#ContType is Temperature
#ContBound = 800

InitialCond {{
	t = {temperature}
{init_condition}
}}
	""".format(tend = target.add["total_time"],tres = target.add["residenceTime"], pressure = target.pressure, init_condition = target.initialCond , temperature = target.temperature, phi = target.phi, simulation = target.simulation, fuel_is = target.fuel_is, g_r = global_reaction,heat = target.add["heat"],T_amb = target.add["T_amb"],isotherm = target.add["isIsotherm"])
			extract = """#!/usr/bin/python
import numpy as np
import os
import pandas as pd
import scipy as sp
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt

flw_target = "{flw_target}"
anchor = {anchor}
limit = {limit}
flw_method = "{flw_method}"
residenceTime = {residenceTime}
anchor_time = {anchor_time}
limit_units = "{limit_units}"
initial_species_conc = {xo}

os.chdir("output")
list_files = os.listdir()

for i in list_files:
	if i.startswith("X"):
		outfile = i
		break
	else:
		continue

#print(outfile)
data = []
read = open(outfile,"r").readlines()
for line in read:
	if "*" not in line:
		data.append(line.strip("\\n"))
#print(data)
file_dump = open("modified.csv","w+")
file_dump.write(("\\n").join(data))
file_dump.close()
file_dump = open("../time_history.csv").write(data)


df = pd.read_csv("modified.csv",sep="\\t")
for i in df:
	if "X-"+flw_target+" " in i:
		profile = df[i]
	if "t" in i:
		time = df[i]/1000
#print(time/1000)
#print(profile)
time_scale = []
for i in time:
	if i<=residenceTime:
		time_scale.append(i)
	elif i == residenceTime:
		time_scale.append(i)
	else:
		continue
time_range = np.asarray(time_scale)
#print(time_range)
def interpolation2D(d, x):
	output = d[0][1] + (x - d[0][0]) * ((d[1][1] - d[0][1])/(d[1][0] - d[0][0]))
	return output

def getNearest(value,array):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return idx
#create a def to populate the array:
def populate(target,array,x):
	new = []
	for i,ele in enumerate(array):
		if i<=x:
			new.append(ele)
			continue
		elif i==x+1:
			new.append(target)
			new.append(ele)
			continue
		else:
			new.append(ele)
			continue
	new = np.asarray(new)
	return new
#create a def to extract the details from the time array
def getTimeStamp(target,array):
	if target in array:
		for i in array:
			if i == target:
				array = list(array)
				x = array.index(i)
				y = len(array[x+1:])
				array = np.asarray(array)
	else:	   
		idx = getNearest(target,array)
		if array[idx] > target:
		#between idx-1 and idx
			x = idx-1
			y = len(array[x:])
			new = populate(target,array,x)
			array = new	
		#can use chebychev interpolation technique
		#start with linear interpolation
		else: 
		#array[idx] < target:
		#between idx and idx+1
			x = idx			
			y = len(array[x:])
			new = populate(target,array,x)
			array = new
	return array,x,y

def timeScaleShift(target,array,x,y):
	#print(x)
	#print(y)
	if target in array:
		for i in array:
			if i == target:
				array = list(array)
				mid = array.index(i)
				#print(mid)
				middle = [array[mid]]
				top = array[mid-x:mid]
				bottom = array[mid+1:mid+y+1]
				scale_shift = np.asarray(top+middle+bottom)
	else:		 
		idx = getNearest(target,array)
		if array[idx] > target:
		#between idx-1 and idx
			x_ = idx-1
			y_ = idx
			new = populate(target,array,x_-1)
			array = new
			scale_shift = timeScaleShift(target,array,x,y)
#		 
#can use chebychev interpolation technique
#		 #start with linear interpolation
		else: 
#		 #array[idx] < target:
#		 #between idx and idx+1
			x_ = idx
			y_ = idx+1
			new = populate(target,array,x_-1)
			#print(new)
			array = new
			scale_shift = timeScaleShift(target,array,x,y)
	#print(len(scale_shift))		#			 array = new
	return scale_shift

def rateDefination(g1,g2):
	t1= g1[0]
	x1 =g1[1]
	t2 = g2[0]
	x2 = g2[1]
	rate = np.abs((x2-x1)/((t1-t2)*1000))
	dt = abs(t2-t1)*1000
	#rate calculated in ppm/ms
	#time calculated in ms
	return rate,dt

def getDataPoints(arr_t,arr_x,X):
	if X in arr_x:
		for i,x in enumerate(arr_x):
			if x == X:
				x_p = x
				t_p = arr_t[i]
	else:
		idx = getNearest(X,arr_x)
	   
		if arr_x[idx] > X:
			x1_ = idx-1
			x2_ = idx
			t1 = arr_t[idx-1]
			t2 = arr_t[idx]
			t_p = interpolation2D(np.array([[arr_x[x1_],t1],[arr_x[x2_],t2]]),X)
			x_p = X  
#		 
#can use chebychev interpolation technique
#		 #start with linear interpolation
		else: 
#		 #array[idx] < target:
#		 #between idx and idx+1
			x1_ = idx
			x2_ = idx+1
			t1 = arr_t[idx]
			t2 = arr_t[idx+1]
			t_p = interpolation2D(np.array([[arr_x[x1_],t1],[arr_x[x2_],t2]]),X)
			x_p = X  
	return (t_p,x_p)   
	
def getRate(xo,t_half,X,t,add):
	if "slope" in add["method"]:	
		if "percentage" in add["unit"]:
			fact = 100
		else:
			fact = 1
		i_y1 = float(add["range_"][0])/fact
		i_y2 = float(add["range_"][1])/fact
		i_anchor = float(add["anchor"])/fact
		i_array_x = X
		i_array_t = t
		ni_array_t,len_x,len_y = getTimeStamp(t_half,i_array_t)
		#print(len(ni_array_t))
		ni_array_x = timeScaleShift(float(xo)*i_anchor,i_array_x,len_x,len_y)
		#print(len(ni_array_x))
		p1 = getDataPoints(i_array_t,ni_array_x,float(xo)*i_y1)
		p2 = getDataPoints(i_array_t,ni_array_x,float(xo)*i_y2)
		rate,dt = rateDefination(p1,p2)
	return ni_array_t,ni_array_x,rate,dt

def load(species,string):
	string = str(string).strip("{{}}")
	st = string.split(",")
	for i in st:
		k = i.split(":")
		if str(species) == k[0].strip('""').strip("''"):
			return_val = float(k[1])
	return return_val
	
add = dict(method = flw_method,range_ = (40,60), anchor = 50,unit = limit_units)
ni_array_t,ni_array_x,rate,dt = getRate(initial_species_conc,anchor_time,profile,time_range,add)

#print(ni_array_t)
#print(ni_array_x)
rate_ = open("rate.csv","w")
stringR = "rate	{{}}	ppm/ms".format(float(rate*1000000))
rate_.write(stringR)
rate_.close()
time_ = open("time.csv","w")
stringT = "time	{{}}	ms".format(dt)
time_.write(stringT)
time_.close()
#fig = plt.figure()
#plt.plot(ni_array_t,ni_array_x,"-")
#plt.show()""".format(flw_target = target.add["flw_species"],anchor = target.add["anchor"],limit = target.add["flw_limits"],flw_method = target.add["flw_method"],residenceTime = target.add["residenceTime"],anchor_time = target.add["anchor_time"],limit_units = target.add["limit_units"],xo = target.species_dict[str(target.add["flw_species"])])
	else:
		raise AssertionError("Unknown type of experiment")	
	
	if target.add["solver"] == "cantera":
		s_convert = """#!/bin/bash
ck2yaml --input=mechanism.mech --thermo=thermo.therm --transport=transport.trans &> out """
		
		s_run = """#!/bin/bash
python3.9 cantera_.py &> solve"""

		extract = """
		"""
	elif target.add["solver"] == "CHEMKIN_PRO":
		s_convert = """#!/bin/bash
ck2yaml --input=mechanism.mech --thermo=thermo.therm --transport=transport.trans &> out """
		s_run = """#!/bin/bash
export NUMEXPR_MAX_THREADS=1
python3.9 {MUQSAC}/soln2ck_2.py --mechanism=mechanism.inp --thermo=thermo.dat {mech} &> soln2ck.out
python3.9 cantera_.py &> solve
		""".format(MUQSAC = opt_dict["Bin"]["bin"] ,Bin=opt_dict["Bin"]["solver_bin"],mech = mech_file,thermo_file = thermo_file_location, trans_file = trans_file_location, fsc = file_specific_command)
		extract = """
		"""
	else:
		s_convert = """#!/bin/bash
python3.9 /home/krithika/Desktop/KineticMechanismOptimization/sc_v2/v2.1/soln2ck.py &>sol2yaml_out
{Bin}/ScanMan -i mechanism.mech -t mechanism.therm -m mechanism_tranport.dat {fsc} -3sr -N 0.05 -E -o mechanism.pre &> ScanMan.log
rm -f ScanMan.log
		""".format(Bin=opt_dict["Bin"]["solver_bin"],thermo_file = thermo_file_location, trans_file = trans_file_location, fsc = file_specific_command) 

			
		s_run = """#!/bin/bash
python3.9 {MUQSAC}/soln2ck_FM.py --mechanism=mechanism.inp --thermo=thermo.dat --transport=transport.dat {mech} &> soln2ck.log
{Bin}/ScanMan -i mechanism.inp -t thermo.dat -m transport.dat -f chemkin -3sr -N 0.05 -E -o mechanism.pre &> ScanMan.log
{Bin}/FlameMan &> Flame.log && python3.9 extract.py &> extract.log
		""".format(MUQSAC = opt_dict["Bin"]["bin"] ,Bin=opt_dict["Bin"]["solver_bin"],mech = mech_file,thermo_file = thermo_file_location, trans_file = trans_file_location, fsc = file_specific_command)
		
	return instring,s_convert,s_run,extract
