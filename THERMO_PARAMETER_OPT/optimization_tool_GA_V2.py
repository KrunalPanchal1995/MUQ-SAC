import numpy as np
from solution import Solution
from scipy.optimize import minimize 
from scipy import optimize as spopt
import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import pygad
import matplotlib as mpl
mpl.use('Agg')
style.use("fivethirtyeight")
import os
from scipy.optimize import rosen, differential_evolution
from scipy.optimize import NonlinearConstraint, Bounds
from scipy.optimize import shgo
from scipy.optimize import BFGS
import pickle
import queue
import threading
import copy

# Global counters for tracking progress and errors
total_evaluations_count = 0
valid_solutions_count = 0
constraint_errors_count = 0

# Lock to ensure thread-safe access to the counters
counter_lock = threading.Lock()



class OptimizationTool(object):
	def __init__(self,
		target_list=None,frequency = None,opt_method = "GA",unsrt = None):
		
		self.target_list = target_list
		self.objective = 0
		self.frequency = frequency
		self.count = 0
		self.opt_method = opt_method
		self.unsrt = unsrt
	
	def obj_func_of_selected_PRS_thermo(self,Z):	
		
		self.count +=1
		string_x = "{self.count},"
		for i in Z:
			string_x+=f"{i},"
		string_x+="\n"
		note = open("guess_values.txt","+a").write(string_x)

		###Algorithm:
		#1] Get the guess values into saperate bins corresponding to species
			# eg: Z_{Ar: Low} = [x1,x2,x3,x4,x5] -> store x5 in self.unsrt[Ar:High].first
			# pass Z_{Ar: Low} to DECODER, get X_{Ar:Low}
		V_of_Z = {}
		count = 0
		for ind,speciesID in enumerate(self.unsrt):
			if self.unsrt[speciesID].a1_star is None or np.all(self.unsrt[speciesID].a1_star) == None:
				species_name = self.unsrt[speciesID].species
				lim = self.unsrt[speciesID].temp_limit
				p_o = self.species_dict[species_name][lim]
				zeta_max = self.zeta_dict[species_name][lim]
				#print(zeta_max)
				
				cov = self.cov_dict[species_name][lim]
				Z_species = Z[count:count+len(p_o)]
				T_low = self.T[count:count+len(p_o)]
				
				Cp_T = []
				for index in range(5):
					t = T_low[index]
					Theta = np.asarray([t/t,t,t**2,t**3,t**4])
					p_ = p_o + Z_species[index]*np.asarray(np.dot(cov,zeta_max)).flatten()
					Cp_T.append(Theta.dot(p_))
				x_species = self.unsrt[speciesID].DECODER(np.asarray(Cp_T),np.asarray(T_low),species_name,self.count,tag="Low")
				count+=len(p_o)
				V_of_Z[speciesID] = x_species
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_star = Z_species
						self.unsrt[spec].a2_star = Z_species[-1]
						self.unsrt[spec].Cp_T_mid_star = Cp_T[-1]
			else:
				species_name = self.unsrt[speciesID].species
				lim = self.unsrt[speciesID].temp_limit
				p_o = self.species_dict[species_name][lim]
				zeta_max = self.zeta_dict[species_name][lim]
				cov = self.cov_dict[species_name][lim]
				Z_species = Z[count:count+len(p_o)-1]
				Cp_T = [self.unsrt[speciesID].Cp_T_mid_star]  ##Start with mid Cp value 1000k
				T_high = self.T[count:count+len(p_o)-1]
				for index in range(4):  ### Start from the next point after 1000 K
					t = T_high[index]
					Theta = np.asarray([t/t,t,t**2,t**3,t**4])
					p_ = p_o + Z_species[index]*np.asarray(np.dot(cov,zeta_max)).flatten()
					Cp_T.append(Theta.dot(p_))
				x_species = self.unsrt[speciesID].DECODER(np.asarray(Cp_T),np.asarray(T_high),species_name,self.count,tag="High")
				count+=len(p_o)-1
				V_of_Z[speciesID] = x_species
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_star = None
						self.unsrt[spec].a2_star = None
						self.unsrt[spec].Cp_T_mid_star = None
		
		string = ""
		V = []
		for spec in self.unsrt:
			V.extend(list(V_of_Z[spec]))
			for k in V_of_Z[spec]:
				string+=f"{k},"
			#x_transformed.extend(temp)
		string+=f"\n"
		V = np.asarray(V)
		#print(V)
		zeta_file = open("zeta_guess_values.txt","+a").write(string)	
		#2] Use X_{Ar:Low} to generate obj. function value from response surface
			#eta_{case_0}= eta([X_{Ar:Low},X_{Ar:High},X_{H2:Low},...])
			#Make sure V_{case_0} = [X_{Ar:Low},X_{Ar:High},X_{H2:Low},...] is in the order of self.unsrt
		"""
		Just for simplicity
		"""
		x = V
		
		obj = 0.0
		rejected_PRS = []
		rejected_PRS_index = []
		target_value = []
		target_stvd = []
		direct_target_value = []
		direct_target_stvd = []
		target_value_2 = []
		case_stvd = []
		case_systematic_error = []
		response_value = []
		response_stvd = []	
		target_weights = []	
		COUNT_Tig = 0
		COUNT_Fls = 0
		COUNT_All = 0	
		COUNT_Flw = 0
		frequency = {}
		diff = []
		diff_2 = []
		diff_3 = {}
		for i,case in enumerate(self.target_list):
			if self.ResponseSurfaces[i].selection == 1:	
				if case.target == "Tig":
					if case.d_set in frequency:
						frequency[case.d_set] += 1
					else:
						frequency[case.d_set] = 1
					COUNT_All +=1
					COUNT_Tig +=1
					#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
					#dataset_weights = (1/COUNT_All)
					
					val = self.ResponseSurfaces[i].evaluate(x)
					#print(val)
					#val,grad = case.evaluateResponse(x)
					f_exp = np.log(case.observed*10)
					#print(f_exp,val)
					#w = 1/(np.log(case.std_dvtn*10))
					w_ = case.std_dvtn/case.observed
					w = 1/w_
					diff.append((val-f_exp)*w)
					#diff.append((val - f_exp)/f_exp)
					diff_2.append(val - f_exp)
					diff_3[case.uniqueID] = (val-f_exp)/f_exp
					response_value.append(val)
					#response_stvd.append(grad)
					#diff.append(val - np.log(case.observed*10))
					target_value.append(np.log(case.observed*10))
					target_value_2.append(np.log(case.observed))
					target_stvd.append(1/(np.log(case.std_dvtn*10)))
					case_stvd.append(np.log(case.std_dvtn*10))
					case_systematic_error.append(abs(np.log(case.observed*10)-val))
					#target_weights.append(dataset_weights)				
				elif case.target == "Fls":
					if case.d_set in frequency:
						frequency[case.d_set] += 1
					else:
						frequency[case.d_set] = 1
					COUNT_All +=1
					COUNT_Fls +=1
					#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
					#dataset_weights = (1/COUNT_All)
					
					val = np.exp(self.ResponseSurfaces[i].evaluate(x))
					response_value.append(val)
					f_exp = case.observed
					w = 1/(case.std_dvtn)
					diff.append((val - f_exp)*w)
					#response_stvd.append(grad)
					target_value.append(case.observed)
					target_value_2.append(case.observed)
					target_stvd.append(1/(case.std_dvtn))
					case_stvd.append(case.std_dvtn)
					case_systematic_error.append(abs(case.observed)-val)
					#target_weights.append(dataset_weights)	
				elif case.target == "Flw":
					if case.d_set in frequency:
						frequency[case.d_set] += 1
					else:
						frequency[case.d_set] = 1
					COUNT_All +=1
					COUNT_Flw +=1
					#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
					#dataset_weights = (1/COUNT_All)
					
					val = self.ResponseSurfaces[i].evaluate(x)
					response_value.append(val)
					f_exp = case.observed
					w = 1/(case.std_dvtn)
					diff.append((val - f_exp)*w)
					#response_stvd.append(grad)
					target_value.append(np.log(case.observed))
					target_value_2.append(np.log(case.observed))
					target_stvd.append(1/(np.log(case.std_dvtn)+abs(np.log(case.observed)-val)))
					case_stvd.append(np.log(case.std_dvtn))
					case_systematic_error.append(abs(np.log(case.observed)-val))
					#target_weights.append(dataset_weights)	
			
		
		diff = np.asarray(diff)
		multiplicating_factors = []
		#multiplicating_factors = np.asarray(target_weights)*np.asarray(target_stvd)
		for i,case in enumerate(self.target_list):
			if self.ResponseSurfaces[i].selection == 1:	
				if case.target == "Tig":
					multiplicating_factors.append(1/COUNT_Tig)
		
				elif case.target == "Fls":
					multiplicating_factors.append(1/COUNT_Fls)
		
		multiplicating_factors= np.asarray(multiplicating_factors)		
		#Giving all datapoints equal weights
		#multiplicating_factor = 1/COUNT_All
		
		for i,dif in enumerate(diff):
			#obj+= multiplicating_factors[i]*(dif)**2
			#obj+= multiplicating_factor*(dif)**2
			#obj+= multiplicating_factors[i]*(dif)**2
			obj+=dif**2		
		
		Diff_3 = open("Dataset_based_obj","+a").write(f"{diff_3}\n")
		get_opt = open("Objective.txt","+a").write(f"{obj}\n")
		string_v = f"{self.count},"
		for i in x:
			string_v+=f"{x},"
		string_v+="\n"
		note = open("guess_values_TRANSFORMED.txt","+a").write(string_v)
		#string_response="#CaseID,Obs,Xvector...\n"
		#for i in range(len(target_value)):
		#	string_response+=f"{target_value[i]},"
		string_response = ""
		for i in range(len(response_value)):
			string_response+=f"{response_value[i]},"
		string_response+="\n"	
		get_target_value = open("response_values.txt","+a").write(string_response)
		#print(obj)
		#if self.count > 3:
		#	raise AssertionError("Stop!")
		#rint(obj)
		return obj
		
	def prepare_z_for_constraints_passing(self,Z):
		V_of_Z = {}
		count = 0
		for ind,speciesID in enumerate(self.unsrt):
			if self.unsrt[speciesID].a1_star is None or np.all(self.unsrt[speciesID].a1_star) == None:
				species_name = self.unsrt[speciesID].species
				lim = self.unsrt[speciesID].temp_limit
				p_o = self.species_dict[species_name][lim]
				zeta_max = self.zeta_dict[species_name][lim]
				#print(zeta_max)
				
				cov = self.cov_dict[species_name][lim]
				Z_species = Z[count:count+len(p_o)]
				T_low = self.T[count:count+len(p_o)]
				
				Cp_T = []
				for index in range(5):
					t = T_low[index]
					Theta = np.asarray([t/t,t,t**2,t**3,t**4])
					p_ = p_o + Z_species[index]*np.asarray(np.dot(cov,zeta_max)).flatten()
					Cp_T.append(Theta.dot(p_))
				x_species = self.unsrt[speciesID].DECODER(np.asarray(Cp_T),np.asarray(T_low),species_name,self.count,tag="Low")
				count+=len(p_o)
				V_of_Z[speciesID] = x_species
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_star = Z_species
						self.unsrt[spec].a2_star = Z_species[-1]
						self.unsrt[spec].Cp_T_mid_star = Cp_T[-1]
			else:
				species_name = self.unsrt[speciesID].species
				lim = self.unsrt[speciesID].temp_limit
				p_o = self.species_dict[species_name][lim]
				zeta_max = self.zeta_dict[species_name][lim]
				cov = self.cov_dict[species_name][lim]
				Z_species = Z[count:count+len(p_o)-1]
				Cp_T = [self.unsrt[speciesID].Cp_T_mid_star]  ##Start with mid Cp value 1000k
				T_high = self.T[count:count+len(p_o)-1]
				for index in range(4):  ### Start from the next point after 1000 K
					t = T_high[index]
					Theta = np.asarray([t/t,t,t**2,t**3,t**4])
					p_ = p_o + Z_species[index]*np.asarray(np.dot(cov,zeta_max)).flatten()
					Cp_T.append(Theta.dot(p_))
				x_species = self.unsrt[speciesID].DECODER(np.asarray(Cp_T),np.asarray(T_high),species_name,self.count,tag="High")
				count+=len(p_o)-1
				V_of_Z[speciesID] = x_species
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_star = None
						self.unsrt[spec].a2_star = None
						self.unsrt[spec].Cp_T_mid_star = None
		
		string = ""
		V = []
		for spec in self.unsrt:
			V.extend(list(V_of_Z[spec]))
			for k in V_of_Z[spec]:
				string+=f"{k},"
			#x_transformed.extend(temp)
		string+=f"\n"
		V = np.asarray(V)		
		return (V)	
	
		
	
	def cons_6_derivative_positive(self, z):
		_species_order_enforced = [
			"AR_Low", "AR_High", "H2_Low", "H2_High", "O_Low", "O_High",
			"O2_Low", "O2_High", "H2O_Low", "H2O_High", "CO_Low", "CO_High",
			"C_Low", "C_High", "HCO_Low", "HCO_High", "OH*_Low", "OH*_High",
			"H_Low", "H_High"
		]
		print("line 342 z value for constraints check \n ", z)
		
		# Fixed temperature points based on your provided ranges.
		T_low_points = np.array([300, 500, 700, 900])
		T_high_points = np.array([1000, 1200, 1400, 1600])

		all_gradients = []
		
		# Ensure the input vector `z` has the correct number of parameters.
		if len(z) != 100:
			raise ValueError("Input vector 'z' must contain exactly 100 parameters.")
		
		start_index = 0
		# Iterate through each species and its corresponding temperature range.
		for speciesID in _species_order_enforced:
			# If the species is 'O_Low' or 'O_High', skip the constraint check.
			if speciesID == "O_Low" or speciesID == "O_High":
				# Simply increment the start_index to move to the next 5 parameters.
				start_index += 5
			else:
				# For all other species, perform the constraint calculation as planned.
				
				# Determine if the species corresponds to the 'Low' or 'High' temperature range.
				temp_limit_key = "Low" if "Low" in speciesID else "High"
				# Extract the base species name (e.g., 'AR' from 'AR_Low').
				species_name = speciesID.split("_")[0]

				# Slice the `z` vector to get the 5 parameters for the current species.
				z_species = z[start_index : start_index + 5]
				
				# Retrieve the nominal parameters (p_o) and covariance matrix (cov) for the species.
				p_o = self.species_dict.get(species_name, {}).get(temp_limit_key)
				cov = self.cov_dict.get(species_name, {}).get(temp_limit_key)
				
				if p_o is None or cov is None:
					raise KeyError(f"Missing data for species '{speciesID}'.")
					
				# Select the correct temperature array based on the species' temperature range.
				if temp_limit_key == "Low":
					temperatures_points = T_low_points
				else:
					temperatures_points = T_high_points
					
				# Calculate the derivative at each of the selected temperature points.
				for T in temperatures_points:
					# This is the derivative of the basis functions with respect to temperature.
					
					theta_derivative = np.array([0, T / T, 2 * T, 3 * T**2, 4 * T**3])
					
					# Perturb the nominal parameters (p_o) using the genetic algorithm parameters (z_species)
					# and the covariance matrix (cov).
					p_modified = p_o + cov @ z_species
					
					# Calculate the derivative by taking the dot product of the modified parameters
					# and the basis function derivatives.
					#derivative = theta_derivative.dot(z)
					derivative = theta_derivative.dot(p_modified)
					
					
					# Store the result in a list.
					all_gradients.append(derivative)

				# Move the start index to the next group of 5 parameters.
				start_index += 5

		# Return the minimum derivative value found across all species and all temperature points.
		#print(np.min(all_gradients))
		print(np.min(all_gradients))
		#raise AssertionError
		return np.min(all_gradients)		
	
iteration_queue = queue.Queue()
objective_queue = queue.Queue()
valid_iterations_count_lock = threading.Lock()
valid_iterations_count = 0
fitness_threshold = 200
file_writing_thread = None
stop_writing_thread = threading.Event()

def file_writer():
    """Thread to write data from queues to files, ensuring no race conditions."""
    global iteration_queue, objective_queue, stop_writing_thread
    with open("zeta_ga_iteration.txt", "a") as iter_file, \
         open("zeta_ga_objective.txt", "a") as obj_file:
        while not stop_writing_thread.is_set() or not iteration_queue.empty():
            try:
                iteration_data = iteration_queue.get(timeout=1)
                objective_data = objective_queue.get(timeout=1)
                iter_file.write(f"{iteration_data}\n")
                obj_file.write(f"{objective_data}\n")
                iteration_queue.task_done()
                objective_queue.task_done()
            except queue.Empty:
                continue

def start_file_writer_thread():
    """Starts the dedicated file-writing thread."""
    global file_writing_thread, stop_writing_thread
    stop_writing_thread.clear()
    file_writing_thread = threading.Thread(target=file_writer)
    file_writing_thread.daemon = True
    file_writing_thread.start()

def stop_file_writer_thread():
    """Signals the file-writing thread to stop and waits for it to finish."""
    global file_writing_thread, stop_writing_thread
    stop_writing_thread.set()
    if file_writing_thread and file_writing_thread.is_alive():
        file_writing_thread.join()

import copy
import numpy as np

# ... other imports
'''
def fitness_func(ga_instance, solution, solution_idx):
    """
    Calculates the fitness for a given solution.
    This function handles both constraint checking and objective calculation.
    """
    global valid_iterations_count_lock, valid_iterations_count

    # Create a deep copy of the unsrt data for this thread
    thread_safe_unsrt = copy.deepcopy(ga_instance.opt_tool.unsrt)
    
    # Create a new, thread-local OptimizationTool instance
    opt_tool = OptimizationTool(unsrt=thread_safe_unsrt)

    # --- Crucial step: Populate the dictionaries on the new instance ---
    opt_tool.species_dict = {}
    opt_tool.zeta_dict = {}
    opt_tool.cov_dict = {}
    opt_tool.T = []
    for species_descriptor in opt_tool.unsrt:
        base_species = species_descriptor.split(":")[0]

        # Populate T list (same logic as run_optimization_with_selected_PRS_thermo)
        lim = opt_tool.unsrt[species_descriptor].temp_limit
        T_temps = opt_tool.unsrt[species_descriptor].temperatures
        
        if lim == "Low":
            opt_tool.T.extend(list(np.linspace(T_temps[0], T_temps[-1], 5)))
        else:
            opt_tool.T.extend(list(np.linspace(T_temps[0], T_temps[-1], 5)[1:]))

        # Populate dictionaries (same logic)
        if base_species not in opt_tool.species_dict:
            opt_tool.species_dict[base_species] = {}
            opt_tool.zeta_dict[base_species] = {}
            opt_tool.cov_dict[base_species] = {}
        
        opt_tool.species_dict[base_species][species_descriptor.split(":")[1]] = opt_tool.unsrt[species_descriptor].nominal[0:5]
        opt_tool.zeta_dict[base_species][species_descriptor.split(":")[1]] = opt_tool.unsrt[species_descriptor].zeta_max.x
        opt_tool.cov_dict[base_species][species_descriptor.split(":")[1]] = opt_tool.unsrt[species_descriptor].cov
    opt_tool.T = np.asarray(opt_tool.T).flatten()
    opt_tool.count = 0
    
    # --- The rest of your existing fitness_func code ---
    try:
        min_derivative_zeta_itr = opt_tool.prepare_z_for_constraints_passing(solution)
        min_derivative = opt_tool.cons_6_derivative_positive(min_derivative_zeta_itr)
        is_valid = min_derivative > 0
    except (KeyError, ValueError) as e:
        print(f"Constraint check failed: {e}. Solution is invalid.")
        is_valid = False

    # ... and so on

    if is_valid:
        # Constraints are met. Calculate fitness from the objective.
        try:
            obj_value = opt_tool.obj_func_of_selected_PRS_thermo(solution)
            
            # PyGAD maximizes fitness, so fitness is 1 / objective value
            # Add a small epsilon to avoid division by zero if obj_value is 0
            fitness_value = 1.0 / (obj_value + 1e-9)
        except Exception as e:
            print(f"Objective function failed: {e}. Giving a very small fitness.")
            raise AssertionError
            fitness_value = 0.000001
            is_valid = False
            
        if is_valid:
            # Enqueue valid iteration data for asynchronous writing
            iteration_queue.put(solution)
            objective_queue.put(obj_value)
            
            with valid_iterations_count_lock:
                valid_iterations_count += 1
                
            return fitness_value
    
    # If constraints are not met or an error occurred, return a very low fitness.
    return 0.000001
'''

def fitness_func(ga_instance, solution, solution_idx):
    """
    Calculates the fitness for a given solution.
    This function handles both constraint checking and objective calculation
    in a thread-safe manner.
    """
    global total_evaluations_count, valid_solutions_count, constraint_errors_count, counter_lock

    # Update total evaluations count
    with counter_lock:
        total_evaluations_count += 1
        # Print progress every 100 iterations
        if total_evaluations_count % 100 == 0:
            print(f"Total evaluated: {total_evaluations_count} | Valid solutions: {valid_solutions_count} | Constraint errors: {constraint_errors_count}")

    # Create a deep copy of the unsrt data for this thread
    thread_safe_unsrt = copy.deepcopy(ga_instance.opt_tool.unsrt)
    
    # Create a new, thread-local OptimizationTool instance
    opt_tool = OptimizationTool(unsrt=thread_safe_unsrt)

    # Populate the dictionaries on the new instance
    opt_tool.species_dict = {}
    opt_tool.zeta_dict = {}
    opt_tool.cov_dict = {}
    opt_tool.T = []
    for species_descriptor in opt_tool.unsrt:
        base_species = species_descriptor.split(":")[0]
        lim = opt_tool.unsrt[species_descriptor].temp_limit
        T_temps = opt_tool.unsrt[species_descriptor].temperatures
        
        if lim == "Low":
            opt_tool.T.extend(list(np.linspace(T_temps[0], T_temps[-1], 5)))
        else:
            opt_tool.T.extend(list(np.linspace(T_temps[0], T_temps[-1], 5)[1:]))

        if base_species not in opt_tool.species_dict:
            opt_tool.species_dict[base_species] = {}
            opt_tool.zeta_dict[base_species] = {}
            opt_tool.cov_dict[base_species] = {}
        
        opt_tool.species_dict[base_species][species_descriptor.split(":")[1]] = opt_tool.unsrt[species_descriptor].nominal[0:5]
        opt_tool.zeta_dict[base_species][species_descriptor.split(":")[1]] = opt_tool.unsrt[species_descriptor].zeta_max.x
        opt_tool.cov_dict[base_species][species_descriptor.split(":")[1]] = opt_tool.unsrt[species_descriptor].cov
    opt_tool.T = np.asarray(opt_tool.T).flatten()
    opt_tool.count = 0
    
    # Check constraints using a try-except block
    is_valid = False
    try:
        min_derivative_zeta_itr = opt_tool.prepare_z_for_constraints_passing(solution)
        print("zeta solutuin is \n ", solution)
        print("after preparing z we have \n ",min_derivative_zeta_itr )
        min_derivative = opt_tool.cons_6_derivative_positive(solution)
        is_valid = min_derivative > 0
    except Exception as e:
        # A non-broadcastable shape error is a type of constraint check failure.
        # This will catch all errors, including the "operands could not be broadcast" error.
        with counter_lock:
            constraint_errors_count += 1
        print(f"Constraint check failed due to error: {e}. Solution is invalid.")
        is_valid = False

    if is_valid:
        # Constraints are met, calculate fitness from the objective
        try:
            obj_value = opt_tool.obj_func_of_selected_PRS_thermo(solution)
            
            # Use a large negative number for objective values close to zero to avoid division by zero
            fitness_value = 1.0 / (obj_value + 1e-12)
            
            # Increment valid solutions counter
            with counter_lock:
                valid_solutions_count += 1
            
            return fitness_value
        except Exception as e:
            # If objective calculation fails, treat as invalid
            print(f"Objective function failed: {e}. Returning low fitness.")
            return -99999999999  # Very low fitness value
    else:
        # Constraints not met, return a very low fitness value
        return -99999999999


def on_generation(ga_instance):
    """Callback function called after each generation."""
    global valid_iterations_count, fitness_threshold
    
    # Check for the valid iteration count termination condition
    with valid_iterations_count_lock:
        if valid_iterations_count >= 1000000:
            print("Termination condition met: 1,000,000 valid iterations reached.")
            ga_instance.run_completed = True
            return

    # Check for the fitness threshold termination condition
    best_solution_fitness = ga_instance.best_solution()[1]
    if best_solution_fitness >= fitness_threshold:
        print(f"Termination condition met: Fitness {best_solution_fitness} >= {fitness_threshold}.")
        ga_instance.run_completed = True
        return

def run_optimization_with_selected_PRS_thermo(unsrt_data, ResponseSurfaces, Input_data):
    """
    The main function to set up and run the genetic algorithm optimization.

    This function initializes the PyGAD instance with parallel processing and 
    the provided objective and constraint functions.
    """
    global valid_iterations_count, iteration_queue, objective_queue

    # Clear previous run data
    valid_iterations_count = 0
    while not iteration_queue.empty():
        iteration_queue.get()
    while not objective_queue.empty():
        objective_queue.get()

    start_file_writer_thread()

    # OptimizationTool instance is required for the fitness function to access the
    # required data (self.unsrt, self.species_dict, etc.)
    opt_tool = OptimizationTool()
    opt_tool.unsrt = unsrt_data
    opt_tool.ResponseSurfaces = ResponseSurfaces
    opt_tool.Input_data = Input_data
    
    # -------------------------------------------------------------------------
    # CORRECTED PART: Build the dictionaries and lists directly from unsrt_data
    # and attach them to the opt_tool instance, as the original script did.
    # -------------------------------------------------------------------------

# -------------------------------------------------------------------------
    # CORRECTED PART: Build the dictionaries and lists directly from unsrt_data
    # and attach them to the opt_tool instance, as the original script did.
    # -------------------------------------------------------------------------

    opt_tool.species_dict = {}
    opt_tool.zeta_dict = {}
    opt_tool.cov_dict = {}
    opt_tool.T = []
    for species_descriptor in opt_tool.unsrt:
        base_species = species_descriptor.split(":")[0]

        # Populate T list
        lim = opt_tool.unsrt[species_descriptor].temp_limit
        T_temps = opt_tool.unsrt[species_descriptor].temperatures
        
        if lim == "Low":
            opt_tool.T.extend(list(np.linspace(T_temps[0], T_temps[-1], 5)))
        else:
            # Skip the common temperature (1000K)
            opt_tool.T.extend(list(np.linspace(T_temps[0], T_temps[-1], 5)[1:]))

        # Populate dictionaries
        if base_species not in opt_tool.species_dict:
            opt_tool.species_dict[base_species] = {}
            opt_tool.zeta_dict[base_species] = {}
            opt_tool.cov_dict[base_species] = {}
        
        opt_tool.species_dict[base_species][species_descriptor.split(":")[1]] = opt_tool.unsrt[species_descriptor].nominal[0:5]
        
        # ADD PRINT STATEMENTS HERE
        #print("running with 100 thread")
        #print(f"Species: {species_descriptor}")
        #print(f"Shape of zeta_max: {opt_tool.unsrt[species_descriptor].zeta_max.x.shape}")
        #print(f"Shape of cov: {opt_tool.unsrt[species_descriptor].cov.shape}")
        
        opt_tool.zeta_dict[base_species][species_descriptor.split(":")[1]] = opt_tool.unsrt[species_descriptor].zeta_max.x
        opt_tool.cov_dict[base_species][species_descriptor.split(":")[1]] = opt_tool.unsrt[species_descriptor].cov

    opt_tool.T = np.asarray(opt_tool.T).flatten()
    opt_tool.count = 0
    # -------------------------------------------------------------------------
    # END OF CORRECTED PART
    # -------------------------------------------------------------------------
    
    # Create the PyGAD instance
    ga_instance = pygad.GA(
        num_generations=2000,
        num_parents_mating=4,
        fitness_func=fitness_func,
        sol_per_pop=200,
        num_genes=100,
        on_generation=on_generation,
        parallel_processing=["thread", 1],
        gene_type=float,
        init_range_low=-1.0,
        init_range_high=1.0,
    )

    # Attach the OptimizationTool instance to the GA instance for easy access
    ga_instance.opt_tool = opt_tool
    ga_instance.opt_tool_data = unsrt_data

    # Run the genetic algorithm
    ga_instance.run()

    # Stop the file-writing thread after the GA completes
    stop_file_writer_thread()

    # Retrieve the best solution found by the GA
    optimal_parameters, solution_fitness, solution_idx = ga_instance.best_solution()
    
    # Now, use the optimal parameters to perform the final decoding process.
    Z_star = optimal_parameters
    V_of_Z = {}
    count = 0
    
    # The `self` from your original snippet is now `opt_tool`
    for ind, speciesID in enumerate(opt_tool.unsrt):
        if opt_tool.unsrt[speciesID].a1_star is None or np.all(opt_tool.unsrt[speciesID].a1_star) == None:
            species_name = speciesID.split(":")[0] # Correctly get species name
            lim = speciesID.split(":")[1] # Correctly get temperature limit
            p_o = opt_tool.species_dict[species_name][lim]
            zeta_max = opt_tool.zeta_dict[species_name][lim]
            cov = opt_tool.cov_dict[species_name][lim]
            Z_species = Z_star[count:count+len(p_o)]
            T_low = opt_tool.T[count:count+len(p_o)]
            Cp_T = []
            for index in range(5):
                t = T_low[index]
                Theta = np.array([t/t, t, t**2, t**3, t**4])
                p_ = p_o + Z_species[index]*np.asarray(np.dot(cov, zeta_max)).flatten()
                Cp_T.append(Theta.dot(p_))
            
            x_species = opt_tool.unsrt[speciesID].DECODER(np.asarray(Cp_T), np.asarray(T_low), species_name, opt_tool.count, tag="Low")
            count += len(p_o)
            V_of_Z[speciesID] = x_species
            for spec in opt_tool.unsrt:
                if opt_tool.unsrt[spec].species == species_name:
                    opt_tool.unsrt[spec].a1_star = Z_species
                    opt_tool.unsrt[spec].a2_star = Z_species[-1]
                    opt_tool.unsrt[spec].Cp_T_mid_star = Cp_T[-1]
        else:
            species_name = speciesID.split(":")[0] # Correctly get species name
            lim = speciesID.split(":")[1] # Correctly get temperature limit
            p_o = opt_tool.species_dict[species_name][lim]
            zeta_max = opt_tool.zeta_dict[species_name][lim]
            cov = opt_tool.cov_dict[species_name][lim]
            Z_species = Z_star[count:count+len(p_o)-1]
            Cp_T = [opt_tool.unsrt[speciesID].Cp_T_mid_star]
            T_high = opt_tool.T[count:count+len(p_o)-1]
            for index in range(4):
                t = T_high[index]
                Theta = np.array([t/t, t, t**2, t**3, t**4])
                p_ = p_o + Z_species[index]*np.asarray(np.dot(cov, zeta_max)).flatten()
                Cp_T.append(Theta.dot(p_))
            x_species = opt_tool.unsrt[speciesID].DECODER(np.asarray(Cp_T), np.asarray(T_high), species_name, opt_tool.count, tag="High")
            count += len(p_o) - 1
            V_of_Z[speciesID] = x_species

    string = ""
    V_star = []
    for spec in opt_tool.unsrt:
        for k in V_of_Z[spec]:
            string += f"{k},"
            V_star.append(k)
    string += f"\n"
    V_star = np.asarray(V_star)

    with open("OPTIMIZED_ZETA.csv", "w") as f:
        f.write(string)

    # -------------------------------------------------------------------------
    # Return the three required outputs.
    # -------------------------------------------------------------------------
    return np.asarray(Z_star), V_star, opt_tool.cov_dict
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




'''
		
	def cons_6_derivative_positive(self, z):
		"""
		This function checks for two constraints in one pass:
		1. Equality of Cp at 1000K for low and high temperature ranges.
		2. Positivity of all derivatives.
		"""
		_species_order_enforced = [
			"AR_Low", "AR_High", "H2_Low", "H2_High", "O_Low", "O_High",
			"O2_Low", "O2_High", "H2O_Low", "H2O_High", "CO_Low", "CO_High",
			"C_Low", "C_High", "HCO_Low", "HCO_High", "OH*_Low", "OH*_High",
			"H_Low", "H_High"
		]
		
		T_low_points = np.array([300, 500, 700, 900])
		T_high_points = np.array([1000, 1200, 1400, 1600])

		all_gradients = []
		start_index = 0

		if len(z) != 100:
			raise ValueError("Input vector 'z' must contain exactly 100 parameters.")
		
		# Iterate through each species and its corresponding temperature range.
		for speciesID in _species_order_enforced:
			species_name = speciesID.split("_")[0]
			temp_limit_key = speciesID.split("_")[1]

			# The slicing is dependent on the temperature range
			if temp_limit_key == "Low":
				z_species_low = z[start_index : start_index + 5]
				p_o_low = self.species_dict.get(species_name, {}).get("Low")
				cov_low = self.cov_dict.get(species_name, {}).get("Low")

				# Check if this species also has a 'High' range for the equality constraint
				if f"{species_name}_High" in _species_order_enforced:
					# --- Equality Constraint Check at 1000K ---
					T_mid = 1000.0
					theta_mid = np.array([T_mid / T_mid, T_mid, T_mid**2, T_mid**3, T_mid**4])
					
					# Get the parameters for the 'High' range
					p_o_high = self.species_dict.get(species_name, {}).get("High")
					cov_high = self.cov_dict.get(species_name, {}).get("High")
					z_species_high = z[start_index + 4: start_index + 4 + 5]

					p_modified_low = p_o_low + cov_low @ z_species_low
					Cp_low_at_1000k = theta_mid.dot(p_modified_low)
					
					p_modified_high = p_o_high + cov_high @ z_species_high
					Cp_high_at_1000k = theta_mid.dot(p_modified_high)

					if not np.isclose(Cp_low_at_1000k, Cp_high_at_1000k, atol=1e-6):
						return -1.0 # Fail fast if equality constraint is not met
			
				# --- Derivative Constraint Check ---
				if speciesID == "O_Low" or speciesID == "O_High":
					# Skip the derivative check for 'O' species
					pass
				else:
					temperatures_points = T_low_points
					for T in temperatures_points:
						theta_derivative = np.array([0, T/T, 2 * T, 3 * T**2, 4 * T**3])
						p_modified = p_o_low + cov_low @ z_species_low
						derivative = theta_derivative.dot(p_modified)
						all_gradients.append(derivative)
				
				start_index += 5

			else: # temp_limit_key == "High"
				z_species_high = z[start_index : start_index + 5]
				p_o_high = self.species_dict.get(species_name, {}).get("High")
				cov_high = self.cov_dict.get(species_name, {}).get("High")

				# --- Derivative Constraint Check ---
				if speciesID == "O_Low" or speciesID == "O_High":
					# Skip the derivative check for 'O' species
					pass
				else:
					temperatures_points = T_high_points
					for T in temperatures_points:
						theta_derivative = np.array([0, T/T, 2 * T, 3 * T**2, 4 * T**3])
						p_modified = p_o_high + cov_high @ z_species_high
						derivative = theta_derivative.dot(p_modified)
						all_gradients.append(derivative)

				start_index += 5

		return np.min(all_gradients)	'''








