import subprocess
import os

def create_start_profile_input(loc_bin,flame_configuration,flame_master_pre_file,available_profile,available_fuel,available_oxidizer,available_bg1,available_bg2,available_bg3,available_global_rxn,target_profile_fuel,target_oxidizer,target_global_rxn,target_fuel_mole_frac,target_ox_mole_frac,target_bg1,target_bg2,target_bg3,target_bg1_x,target_bg2_x,target_bg3_x,target_pressure,target_temperature):
	instring = ''
	
	instring = '''

#########
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

#######################
# Boundary conditions #
#######################

Boundary:
 
 fuel: {key: %(key_f)s,From: %(from_f)s, To: %(to_f)s,ToConc: %(to_concF)s}
 
 oxidizer: {key: %(key_ox)s, From: %(from_ox)s, To: %(to_ox)s, ToConc : %(to_concOx)s}
 
 bathGas: {key: %(key_bg)s, From: %(from_bg)s, To:  %(to_bg)s,ToConc:%(to_concBg)s}
 
 globalReaction: {FromRxn: %(from_rxn)s, ToRxn: %(to_rxn)s}
 
 pressure: %(to_Pa)s 
 
 temperature: %(to_K)s
 
 units: {pressure: %(unit_p)s, temperature: %(unit_t)s}
	
'''% {"flame_config":flame_configuration,"pre_file":flame_master_pre_file,"init_profile":available_profile,"fuel_a":available_fuel,"oxidizer":available_oxidizer,"bg1_a":available_bg1,"bg2_a":available_bg2,"bg3_a":available_bg3,"globle_rxn_a":available_global_rxn,"fuel_t":target_profile_fuel,"oxidizer_t":target_oxidizer,"globle_rxn_t":target_global_rxn,"fuel_t_x":target_fuel_mole_frac,"oxidizer_t_x":target_ox_mole_frac,"bg1_t":target_bg1,"bg2_t":target_bg2,"bg3_t":target_bg3,"bg1_t_x":target_bg1_x,"bg2_t_x":target_bg2_x,"bg3_t_x":target_bg3_x,"P":target_pressure,"T":target_temperature}

	infile = open("profile_generator.input",'w')
	infile.write(instring)
	infile.close()
	
	run_file = open("run",'w')
	s = """#!/bin/bash
python %(profile_generator_src)s profile_generator.input &> profile.log
"""% {"profile_generator_src":loc_bin}
#&> Flame.log ; gnome-terminal && tail -f flame.log
	run_file.write(s)
	run_file.close()

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
		instring = """
############
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

WriteEverySolution = TRUE
PrintMolarFractions is TRUE

OutputPath is ./Output
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

Unburnt Side {
	dirichlet {{
		t = {temperature}
		{dirichlet}
	}}
}

{bounds}
""".format(time_flag=inputs["TimeFlag"],init_grid=inputs["initialgridpoints"],max_grid=inputs["maxgridpoints"],s_p_loc=inputs["StartProfilesFile"],p_f=inputs["MechanismFile"],g_r=globeRxn,fuel_is=fuel["fuelIs"],oxidizer=oxidizer["oxidizerIs"],flame=inputs["Flame"],pressure=P,radiation_tag=inputs["ComputeWithRadiation"],thermodiffusion_tag=inputs["Thermodiffusion"],etp=inputs["ExpTempFile"],flow_rate=inputs["flowRate"],temperature=T,dirichlet=controls["dirichlet"],bounds=controls["conc_bounds"])

	infile = open("FlameMaster.input",'w')
	infile.write(instring)
	infile.close()
	
	run_file = open("run",'w')
	s = """#!/bin/bash
/home/krunal/FlameMaster/Bin/bin/FlameMan 
"""
#&> Flame.log ; gnome-terminal && tail -f flame.log
	run_file.write(s)
	run_file.close()
	subprocess.call(["chmod","+x",'run'])

def create_input_file(opt_dict, target,fuel, g_reaction, thermo_file_location, trans_file_location,startProfile_location,file_specific_command):
	
	global_reaction = ''
	instring = ''
	inputs = target.inputs
	fuel = target.fuel
	P = target.pressure
	#print(target.target)
	for i in g_reaction:
		global_reaction += i
		
	# Defining the default variables with solver specific keys:
	# In later version this moduel is to be shifted in combustion_target_class
	# creating an input dictionary, lateron the dictionary will be part of target object
	
	default = {}
	cantera_def = {}
	cantera_def["saveAll"] = '#' 
	cantera_def["ign_species"] = 'OH'
	cantera_def["ign_cond"] = 'max'
	cantera_def["ign_specific_cond"] = "None;None"
	FlameMan_def = {}
	

		
	if target.input_file != None:
		instring = open(target.input_file,'r').read()
	
	elif "Tig" in target.target:
		if target.ignition_type == "reflected":
			if target.add["solver"] == "cantera":
				instring = """from __future__ import division
from __future__ import print_function
import pandas as pd
import numpy as np
import time
import cantera as ct

gas = ct.Solution('mech.cti')

reactorTemperature = {temperature} #Kelvin
reactorPressure = {pressure}

gas.TPX = reactorTemperature, reactorPressure,{species_conc}
r = ct.Reactor(contents=gas)
reactorNetwork = ct.ReactorNet([r])

stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)

def ignitionDelay(df, species,cond="max",specified_value ="None;None"):
	if cond == "max":
		tau = df[species].idxmax()
	elif cond == "onset":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(df[species].to_numpy())
		slope = conc/time
		index = int(np.diff(slope).argmax())
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
				index = np.where(np.abs(df[species].to_numpy()-np.ones(len(df[species].to_numpy()))*target)<0.15*target)[0]
				tau = df.index.to_numpy()[index[0]]
			else:
				tau = df[species].eq(target).idxmax()
	return tau

t0 = time.time()
estimatedIgnitionDelayTime = 1
t = 0

counter = 1;
while(t < estimatedIgnitionDelayTime):
    t = reactorNetwork.step()
    #print(t)
    if (counter%1 == 0):
        # We will save only every 1st value (20th value). Otherwise, this takes too long
        # Note that the species concentrations are mass fractions
        timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
    counter+=1

tau = ignitionDelay(timeHistory,{delay_def},{delay_cond},{specific_cond})
# Toc
t1 = time.time()
tau_file = open("tau.out",'w')
tau_file.write("tau 	{{}} us".format(tau*(1000000)))
tau_file.close()
{saveAll}timeHistory.to_csv("time_history.csv")
""".format(temperature = target.temperature,pressure=target.pressure,species_conc = target.species_dict,delay_def = target.ign_def,delay_cond = target.ign_cond,specific_cond = target.ign_specific_cond,saveAll=target.saveAll)
			else:
				instring = """############
	# Numerics #
	############

	RelTol = 1.0e-10
	AbsTol = 1.0e-12

	TStart = 0.0
	TEnd = 1.0e0

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
	#WriteEverySolution is TRUE
	#PrintMolarFractions is TRUE

	OutputPath is ./output
	NOutputs = 50

	#############
	# Chemistry #
	#############

	MechanismFile is mechanism.pre
	globalReaction is {g_r};

	fuel is {f}
	oxidizer is O2

	#########
	# Flame #
	#########

	Flame is {simulation} Homo Reactor
	#Flame is Isobar Homo Reactor
	#ExactBackward is TRUE

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
		X->{f} = {fuel_x}
		X->{ox} = {ox_x}
	  {bg1_tag}	X->{bg1} = {bg1_x}
	  {bg2_tag}	x->{bg2} = {bg2_x}
	  {bg3_tag}	x->{bg3} = {bg3_x}
		
	}}
	""".format(pressure = target.pressure, fuel_x = target.fuel_x,ox = target.oxidizer, ox_x = target.oxidizer_x, bg1_tag = target.bg1_tag, bg1 = target.bath_gas1, bg1_x = target.bath_gas1_x, bg2_tag = target.bg2_tag, bg2 = target.bath_gas2, bg2_x = target.bath_gas2_x, bg3_tag = target.bg3_tag, bg3 = target.bath_gas3, bg3_x = target.bath_gas3_x, temperature = target.temperature, phi = target.phi, simulation = target.simulation, f = fuel, g_r = global_reaction)
		else:
			raise AssertionError("Invalid ignition mode specified for ignition delay simulations")
	elif "Fls" in target.target:
		if target.add["solver"] == cantera:
			instring = """
				
				"""
				
		else:
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
#"time_flag":inputs["TimeFlag"],"init_grid":inputs["initialgridpoints"],"max_grid":inputs["maxgridpoints"],"s_p_loc":inputs["StartProfilesFile"],"p_f":inputs["MechanismFile"],"g_r":globeRxn,"fuel_is":fuel["fuelIs"],"oxidizer":oxidizer["oxidizerIs"],"flame":inputs["Flame"],"phi":phi,"pressure":P,"radiation_tag":inputs["ComputeWithRadiation"],"thermodiffusion_tag":inputs["Thermodiffusion"],"temperature":T,"dirichlet":controls["dirichlet"],"conc_bounds":controls["conc_bounds"])

	elif "Flf" in target.target:
		if target.add["solver"] == cantera:
				
			instring = """
				
				"""
					
		else:
		
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

WriteEverySolution = TRUE
PrintMolarFractions is TRUE

OutputPath is ./Output
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

Unburnt Side {
	dirichlet {{
		t = {temperature}
		{dirichlet}
	}}
}

{bounds}
""".format(time_flag=inputs["TimeFlag"],init_grid=inputs["initialgridpoints"],max_grid=inputs["maxgridpoints"],s_p_loc=inputs["StartProfilesFile"],p_f=inputs["MechanismFile"],g_r=globeRxn,fuel_is=fuel["fuelIs"],oxidizer=oxidizer["oxidizerIs"],flame=inputs["Flame"],pressure=P,radiation_tag=inputs["ComputeWithRadiation"],thermodiffusion_tag=inputs["Thermodiffusion"],etp=inputs["ExpTempFile"],flow_rate=inputs["flowRate"],temperature=T,dirichlet=controls["dirichlet"],bounds=controls["conc_bounds"])
	
	elif "Flw" in target.target:
	#No startProfile needed
		if target.reactor == "JSR" or target.reactor =="PSR" or target.reactor =="FlowReactor":
			if target.add["solver"] == cantera:
				
				instring = """
				
				"""
				
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
	NOutputs = 10000

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
	
	Isothermal is {isotherm}
	
	HeatTransCoeff is {heat}
	
	AmbientTemp is (T_amb)
	
	
	ExactBackward is TRUE

	{phi}

	Pressure = {pressure}


	#######################
	# Boundary conditions #
	#######################

	#ContInc = -25
	#ContType is Temperature
	#ContBound = 800

	InitialCond {{
		t = {temperature}
		{init_cond}
	}}
	{bounds}
	""".format(tend = float(target.Tend/1000),tres = target.additional["TRes"],pressure = target.pressure, temperature = target.temperature, phi = target.phi, simulation = target.simulation, isotherm = target.additional["Isothermal"], heat = target.additional["HeatTransCoeff"], T_amb = target.additional["AmbientTemp"], fuel_is = fuel, g_r = global_reaction, init_cond=controls["init_cond"],bounds=controls["conc_bounds"])
	else:
		raise AssertionError("Unknown type of experiment")	
	
	infile = open("FlameMaster.input",'w')
	infile.write(instring)
	infile.close()
	
	run_file = open("run",'w')
	s = """#!/bin/bash
/home/krunal/FlameMaster/Bin/bin/ScanMan -i mechanism.mech -t {thermo_file} -m {trans_file} {fsc} -o mechanism.pre -3sr &> ScanMan.log
/home/krunal/FlameMaster/Bin/bin/FlameMan &> Flame.log
""".format(thermo_file = thermo_file_location, trans_file = trans_file_location, fsc = file_specific_command) 
	
	run_file.write(s)
	run_file.close()
	#print(os.getcwd())
	subprocess.call(["chmod","+x",'run'])
