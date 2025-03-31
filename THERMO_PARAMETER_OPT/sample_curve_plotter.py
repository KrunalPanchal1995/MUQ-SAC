import os
import numpy as np
import matplotlib.pyplot as plt

def sample_plot(unsrt_data, design_matrix):
	if not os.path.exists("sample_plots"):
		os.mkdir("sample_plots")
	sample_plots_dir = "sample_plots"
	# Access the Lcp and NominalParam attributes from the unsrt data
	Lcp_low = unsrt_data["H2:Low"].cov
	eta_low = unsrt_data["H2:Low"].zeta_max.x	# Access the Lcp attribute
	Lcp_High = unsrt_data["H2:High"].cov
	eta_High = unsrt_data["H2:High"].zeta_max.x
	#print("Lcp", Lcp)
	
	NominalParam_low = np.asarray(unsrt_data["H2:Low"].NominalParams).flatten()  # Convert to array and reshape to 5x1 if needed
	NominalParam_High = np.asarray(unsrt_data["H2:High"].NominalParams).flatten() 
	# Print NominalParam to check its value
	#print("NominalParam:\n", NominalParam)
	#print(Lcp,eta,NominalParam)
	# Define the temperature range
	#T = np.linspace(300, 1000, 20)
	T_low = np.linspace(300,1000,100)  # Temperature range from 298 K to 1000 K
	T_High = np.linspace(1000,5000,100)
	Tscale = 1
	# Calculate theta = [T/T, T, T^2, T^3, T^4] for each temperature
	theta_low = np.array([T_low / T_low, T_low/Tscale, (T_low/Tscale)**2, (T_low/Tscale)**3, (T_low/Tscale)**4]).T  # Transposed to get a 100x5 matrix
	theta_High = np.array([T_High / T_High, T_High/Tscale, (T_High/Tscale)**2, (T_High/Tscale)**3, (T_High/Tscale)**4]).T
	#print(T,theta)
	Cp_max_low = NominalParam_low + Lcp_low.dot( eta_low)
	Cp_min_low = NominalParam_low - Lcp_low.dot(eta_low)
	
	Cp_max_High = NominalParam_High + Lcp_High.dot( eta_High)
	Cp_min_High = NominalParam_High - Lcp_High.dot(eta_High)
	
	
	# Calculate NominalParam * Theta for each temperature (resulting in a 100x1 vector)
	nominal_param_theta_low = theta_low.dot(NominalParam_low)
	nominal_param_theta_High = theta_High.dot(NominalParam_High)
	
	#print(nominal_param_theta)
	cp_param_theta_max_low = 8.314*np.dot(theta_low, Cp_max_low)
	cp_param_theta_min_low = 8.314*np.dot(theta_low, Cp_min_low)
	cp_param_theta_max_High = 8.314*np.dot(theta_High, Cp_max_High)
	cp_param_theta_min_High = 8.314*np.dot(theta_High, Cp_min_High)
	fig = plt.figure(figsize=(10, 6))
	ax = fig.add_subplot()
	# Iterate over each row of the design_matrix to create and save individual plots
	for i, row in enumerate(design_matrix[2052:2102]):
		# Extract the first 5 elements to form zeta (5x1 vector)
		#print(row)
		zeta_low = np.asarray(row[10:15]).reshape(-1, 1)  # Ensure zeta is a 5x1 vector
		zeta_High = np.asarray(row[15:20]).reshape(-1, 1)
		# Compute the result for each temperature
		result_low = 8.314*(np.dot(theta_low, (Lcp_low @ zeta_low)).flatten() + nominal_param_theta_low.flatten())  # Apply matrix operations 
		result_High = 8.314*(np.dot(theta_High, (Lcp_High @ zeta_High)).flatten() + nominal_param_theta_High.flatten())
		#print(result)
		# Plotting the curve
		#print(result_low/8.314)
		#print(result_High/8.314)
		plt.plot(T_low, result_low, color='red',linewidth=0.4)
		plt.plot(T_High, result_High , color = 'green' , linewidth =0.4)
		
	

		# Save the plot to the 'sample_plots' directory
	string_low = ""
	for index,i in enumerate(result_low):
		string_low+=f"{T_low[index]},{result_low[index]}\n"
	string_high = ""
	for index,i in enumerate(result_low):
		string_high+=f"{T_High[index]},{result_High[index]}\n"
	
	f = open("low_T_cp_H2.csv","w").write(string_low)
	g = open("high_T_cp_H2.csv","w").write(string_high)
	plt.plot(T_low, result_low, label='Sample Curve(T<1000K)', color='red',linewidth=0.6)
	plt.plot(T_low, 8.314*nominal_param_theta_low, label='Nominal', color='blue')
	plt.plot(T_High, result_High, label='Sample Curve (T>1000K)', color='green',linewidth=0.6)
	plt.plot(T_High, 8.314*nominal_param_theta_High, label='Nominal', color='blue')
	plt.plot(T_low, cp_param_theta_max_low,'k--')
	plt.plot(T_low, cp_param_theta_min_low,'k--')# ,label='Nominal', color='blue')
	plt.plot(T_High, cp_param_theta_max_High,'k--' ,label='Unsrt. Limits')
	plt.plot(T_High, cp_param_theta_min_High,'k--')# ,label='Nominal', color='blue')
	#plt.title(f'ALL_Zeta Samples for Design Matrix')
	plt.ylim(0.92*min(cp_param_theta_min_low),1.15*max(cp_param_theta_max_High))
	for spine in ax.spines.values():
	    spine.set_edgecolor('black')
	    spine.set_linewidth(1) 
	plt.xlabel('Temperature (K)',fontsize=25)
	plt.ylabel(r'Heat Capacity (Cp, $J mol^{-1} K^{-1}$))',fontsize=25)
	plt.legend(loc="best",fontsize=15)
	ax.grid(False)
	#plt.grid(True)
	plt.savefig(f'sample_plots/H2_Higher_50.pdf',bbox_inches="tight")
	plt.close()  # Close the plot to free up memory

	#print("Plots have been saved to the 'sample_plots' directory.")
	print(f"Sample plots will be saved to: {os.path.abspath(sample_plots_dir)}")

