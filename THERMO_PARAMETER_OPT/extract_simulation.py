'''
# use this code to strore all the entries of simulation from the opt folder
import os
import csv

# Path to the main Opt directory
base_dir = '/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/Opt'
num_cases = 90
num_samples = 36050  # from 0 to 36049

# Output CSV file
output_csv = 'simulation_results.csv'

# Collect all rows (each row corresponds to a case)
all_data = []

for case_num in range(num_cases):
    case_folder = os.path.join(base_dir, f'case-{case_num}')
    row = []
    for sample_num in range(num_samples):
        tau_file = os.path.join(case_folder, str(sample_num), 'output', 'tau.out')
        try:
            with open(tau_file, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    tau_value = lines[1].strip().split()[1]  # Extract second value in second line
                    row.append(tau_value)
                else:
                    row.append('NA')  # Handle missing/invalid content
        except FileNotFoundError:
            row.append('NA')  # Handle missing files
    all_data.append(row)

# Write to CSV
with open(output_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(all_data)

print(f"Extraction complete. CSV saved to: {output_csv}")
'''

'''
import pandas as pd

# Replace with your file path
file_path1 = '/data2/RANA/FINAL_RUN_7nc/DesignMatrix.csv'
file_path2 = '/data2/RANA/FINAL_RUN_7nc/simulation_results.csv'

# Load the CSV file
df = pd.read_csv(file_path1)
df = pd.read_csv(file_path2)


# Print the shape
print(f"Shape of '{file_path1}': {df.shape[0]} rows × {df.shape[1]} columns")
print(f"Shape of '{file_path2}': {df.shape[0]} rows × {df.shape[1]} columns")
'''


# use this to do transport of the above created simulation_result.csv file
import pandas as pd

# Input and output file paths
input_csv = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/design_matrix_zeta_conversion/simulation_results.csv"
output_csv = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/design_matrix_zeta_conversion/simulation_results_transpose.csv"

# Read the CSV file
df = pd.read_csv(input_csv, header=None)

# Transpose the DataFrame
df_transposed = df.T

# Save the transposed DataFrame
df_transposed.to_csv(output_csv, index=False, header=False)

print(f"Transposed CSV saved to '{output_csv}' with shape {df_transposed.shape}")

