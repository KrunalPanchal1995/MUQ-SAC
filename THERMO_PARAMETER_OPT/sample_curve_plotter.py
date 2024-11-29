import os
import numpy as np
import matplotlib.pyplot as plt

def sample_plot(unsrt_data, design_matrix):
	if not os.path.exists("sample_plots"):
	    os.mkdir("sample_plots")
	sample_plots_dir = "sample_plots"
	# Access the Lcp and NominalParam attributes from the unsrt data
	Lcp = unsrt_data["AR:Low"].Lcp  # Access the Lcp attribute
	print("Lcp", Lcp)
	NominalParam = np.asarray(unsrt_data["AR:Low"].NominalParams).reshape(-1, 1)  # Convert to array and reshape to 5x1 if needed
	# Print NominalParam to check its value
	print("NominalParam:\n", NominalParam)

	# Define the temperature range
	T = np.linspace(300, 1000, 100)  # Temperature range from 298 K to 1000 K

	# Calculate theta = [T/T, T, T^2, T^3, T^4] for each temperature
	theta = np.array([T / T, T, T**2, T**3, T**4]).T  # Transposed to get a 100x5 matrix

	# Calculate NominalParam * Theta for each temperature (resulting in a 100x1 vector)
	nominal_param_theta = np.dot(theta, NominalParam)

	# Iterate over each row of the design_matrix to create and save individual plots
	for i, row in enumerate(design_matrix):
	    # Extract the first 5 elements to form zeta (5x1 vector)
	    zeta = np.asarray(row[:5]).reshape(-1, 1)  # Ensure zeta is a 5x1 vector
	    
	    # Compute the result for each temperature
	    result = np.dot(theta, (Lcp @ zeta) + NominalParam.flatten())  # Apply matrix operations to match dimensions
	    print(result)
	    # Plotting the curves
	    plt.figure(figsize=(10, 6))
	    plt.plot(T, result, label='Sample Curve', color='blue')
	    plt.plot(T, nominal_param_theta, label='Nominal', color='red')

	    # Adding plot details
	    plt.title(f'Zeta Samples for Design Matrix Row {i}')
	    plt.xlabel('Temperature (K)')
	    plt.ylabel('Values')
	    plt.legend()
	    plt.grid(True)

	    # Save the plot to the 'sample_plots' directory
	    plt.savefig(f'sample_plots/zeta_samples_row_{i}.png')
	    plt.close()  # Close the plot to free up memory

	#print("Plots have been saved to the 'sample_plots' directory.")
	print(f"Sample plots will be saved to: {os.path.abspath(sample_plots_dir)}")

