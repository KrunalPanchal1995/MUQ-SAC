
import pickle
import os

file_path = '/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/unsrt.pkl'
output_file_path = 'unsrt_data.txt'

if not os.path.exists(file_path):
    print(f"Error: The file '{file_path}' was not found.")
else:
    try:
        with open(file_path, 'rb') as f:
            data = pickle.load(f)

        with open(output_file_path, 'w') as out_f:
            out_f.write("Extracting data for each species:\n")
            out_f.write("-" * 40 + "\n")

            for species_name, species_obj in data.items():
                nominal_value = getattr(species_obj, 'nominal', 'N/A')
                covariance_matrix = getattr(species_obj, 'cov', 'N/A')
                zeta_max = getattr(species_obj, 'zeta_max', 'N/A')

                out_f.write(f"Species: {species_name}\n")
                out_f.write(f"  - Nominal Value: {nominal_value}\n")
                out_f.write(f"  - Covariance Matrix:\n{covariance_matrix}\n")
                out_f.write(f"  - Zeta Max: {zeta_max}\n")
                out_f.write("-" * 40 + "\n")

        print(f"Extraction complete! Data has been saved to '{output_file_path}'.")

    except pickle.UnpicklingError as e:
        print(f"Error: Failed to unpickle the file. Details: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
