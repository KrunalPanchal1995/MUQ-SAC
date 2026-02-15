'''
# use this code to plot cp curvbes for species using 1 given mechanism... each species has one plot
import numpy as np
import matplotlib.pyplot as plt
from ruamel.yaml import YAML

def nasa7_cp(T, coeffs):
    """
    Compute Cp/R from NASA7 coefficients at temperature T
    Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    """
    a1, a2, a3, a4, a5, a6, a7 = coeffs
    return a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4

def load_species_data(yaml_path):
    yaml = YAML(typ='safe')
    with open(yaml_path, 'r') as f:
        data = yaml.load(f)
    return data['species']



import os

def plot_cp_for_species(yaml_path, species_list):
    species_data = load_species_data(yaml_path)

    # make output folder if not exists
    outdir = "opt_cp_curves"
    os.makedirs(outdir, exist_ok=True)

    Tmin = 200
    Tmax = 3000
    T = np.linspace(Tmin, Tmax, 300)

    for sp in species_list:
        found = False
        for entry in species_data:
            if entry['name'].lower() == sp.lower():
                found = True
                trange = entry['thermo']['temperature-ranges']
                coeffs_low = entry['thermo']['data'][0]
                coeffs_high = entry['thermo']['data'][1]

                cp_vals = []
                for temp in T:
                    if temp <= trange[1]:
                        cp = nasa7_cp(temp, coeffs_low)
                    else:
                        cp = nasa7_cp(temp, coeffs_high)
                    cp_vals.append(cp * 8.314)  # convert Cp/R to Cp [J/mol-K]

                # plot & save for each species
                plt.figure()
                plt.plot(T, cp_vals, label=sp)
                plt.xlabel("Temperature [K]")
                plt.ylabel("Cp [J/mol-K]")
                plt.title(f"Cp vs T for {sp}")
                plt.legend()
                plt.grid(True)
                plt.tight_layout()
                plt.savefig(os.path.join(outdir, f"{sp}_Cp.png"), dpi=300)
                plt.close()

        if not found:
            print(f"Species {sp} not found in file.")



if __name__ == "__main__":
    # <<< Put your file path and species names here >>>
    yaml_file = "/data/TEST-THERMO-sens/FFCM1.0/FFCM1.yaml"
    species_to_plot = ["AR", "H2", "O2", "H2O" , "C" , "CO", "OH*", "O" , "HCO" , "H"]  # <-- manually change this list
    plot_cp_for_species(yaml_file, species_to_plot)
'''

####################################################################################
'''
# use this to generate compare plots of cp for origional vs optimized...one figure for both origional and optimized species cp curve

import numpy as np
import matplotlib.pyplot as plt
from ruamel.yaml import YAML
import os

def nasa7_cp(T, coeffs):
    """
    Compute Cp/R from NASA7 coefficients at temperature T
    Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    """
    a1, a2, a3, a4, a5, a6, a7 = coeffs
    return a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4

def load_species_data(yaml_path):
    yaml = YAML(typ='safe')
    with open(yaml_path, 'r') as f:
        data = yaml.load(f)
    return data['species']

def get_cp_curve(species_data, sp, T):
    """Return Cp curve for a given species over T range"""
    for entry in species_data:
        if entry['name'].lower() == sp.lower():
            trange = entry['thermo']['temperature-ranges']
            coeffs_low = entry['thermo']['data'][0]
            coeffs_high = entry['thermo']['data'][1]

            cp_vals = []
            for temp in T:
                if temp <= trange[1]:
                    cp = nasa7_cp(temp, coeffs_low)
                else:
                    cp = nasa7_cp(temp, coeffs_high)
                cp_vals.append(cp * 8.314)  # convert Cp/R to Cp [J/mol-K]
            return cp_vals
    return None

def plot_cp_comparison(yaml_orig, yaml_opt, species_list):
    species_orig = load_species_data(yaml_orig)
    species_opt = load_species_data(yaml_opt)

    outdir = "opt_cp_curves_compare"
    os.makedirs(outdir, exist_ok=True)

    Tmin = 200
    Tmax = 3000
    T = np.linspace(Tmin, Tmax, 300)

    for sp in species_list:
        cp_orig = get_cp_curve(species_orig, sp, T)
        cp_opt = get_cp_curve(species_opt, sp, T)

        if cp_orig is None:
            print(f"Species {sp} not found in original file.")
            continue
        if cp_opt is None:
            print(f"Species {sp} not found in optimized file.")
            continue

        # plot
        plt.figure()
        plt.plot(T, cp_orig, 'r-', label=f"{sp} (Original)")
        plt.plot(T, cp_opt, 'g--', label=f"{sp} (Optimized)")
        plt.xlabel("Temperature [K]")
        plt.ylabel("Cp [J/mol-K]")
        plt.title(f"Cp vs T for {sp}")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{sp}_Cp_compare.png"), dpi=300)
        plt.close()


if __name__ == "__main__":
    yaml_file_orig = "/data/TEST-THERMO-sens/FFCM1.0/FFCM1.yaml"
    yaml_file_opt  = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/GA_ANALYSIS_with_opt_zeta/opt_itr_255885/optimized_results/OPTIMIZED_MECH/new_mech.yaml"

    species_to_plot = ["AR", "H2", "O2", "H2O", "C", "CO", "OH*", "O", "HCO", "H"]
    plot_cp_comparison(yaml_file_orig, yaml_file_opt, species_to_plot)
'''
#########################################################################################################################
'''
# usr this code to plot origional and optimized cp curves along with the uncertainty band (grid wise)

import numpy as np
import matplotlib.pyplot as plt
from ruamel.yaml import YAML
import os
import pickle

def nasa7_cp(T, coeffs):
    """
    Compute Cp/R from NASA7 coefficients at temperature T
    Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T**4
    """
    coeffs = np.array(coeffs)
    a1, a2, a3, a4, a5, a6, a7 = coeffs
    return a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4

def load_species_data(yaml_path):
    yaml = YAML(typ='safe')
    with open(yaml_path, 'r') as f:
        data = yaml.load(f)
    return data['species']

def get_cp_curve(species_data, sp, T):
    """Return Cp curve for a given species over T range"""
    for entry in species_data:
        if entry['name'].lower() == sp.lower():
            trange = entry['thermo']['temperature-ranges']
            coeffs_low = np.array(entry['thermo']['data'][0])
            coeffs_high = np.array(entry['thermo']['data'][1])

            cp_vals = []
            for temp in T:
                if temp <= trange[1]:
                    cp = nasa7_cp(temp, coeffs_low)
                else:
                    cp = nasa7_cp(temp, coeffs_high)
                cp_vals.append(cp * 8.314)
            return np.array(cp_vals)
    return None

def load_uncertainty_data(pkl_path):
    """Load uncertainty data from the pkl file"""
    if not os.path.exists(pkl_path):
        print(f"Error: The file '{pkl_path}' was not found.")
        return None
    try:
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
        return data
    except pickle.UnpicklingError as e:
        print(f"Error: Failed to unpickle the file. Details: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while loading the pkl file: {e}")
        return None

def get_uncertainty_band_curves(species_orig_data, uncertainty_data, sp, T):
    """
    Calculates the upper and lower Cp curves based on uncertainty data.
    """
    orig_entry = None
    for entry in species_orig_data:
        if entry['name'].lower() == sp.lower():
            orig_entry = entry
            break
    
    if not orig_entry:
        print(f"Species {sp} not found in original data.")
        return None, None

    trange = orig_entry['thermo']['temperature-ranges']
    coeffs_low_orig = np.array(orig_entry['thermo']['data'][0])
    coeffs_high_orig = np.array(orig_entry['thermo']['data'][1])

    try:
        nominal_obj_low = uncertainty_data[f"{sp}:Low"]
        nominal_obj_high = uncertainty_data[f"{sp}:High"]
        
        cov_low = getattr(nominal_obj_low, 'cov', None)
        zeta_max_obj_low = getattr(nominal_obj_low, 'zeta_max', None)
        zeta_max_vec_low = getattr(zeta_max_obj_low, 'x', None)
        
        cov_high = getattr(nominal_obj_high, 'cov', None)
        zeta_max_obj_high = getattr(nominal_obj_high, 'zeta_max', None)
        zeta_max_vec_high = getattr(zeta_max_obj_high, 'x', None)
    except (KeyError, AttributeError):
        print(f"Uncertainty data for {sp} is incomplete or not in expected format.")
        return None, None

    delta_low_coeffs = cov_low @ zeta_max_vec_low
    delta_high_coeffs = cov_high @ zeta_max_vec_high

    upper_coeffs_low = coeffs_low_orig.copy()
    lower_coeffs_low = coeffs_low_orig.copy()
    upper_coeffs_high = coeffs_high_orig.copy()
    lower_coeffs_high = coeffs_high_orig.copy()

    upper_coeffs_low[0:5] += delta_low_coeffs
    lower_coeffs_low[0:5] -= delta_low_coeffs
    upper_coeffs_high[0:5] += delta_high_coeffs
    lower_coeffs_high[0:5] -= delta_high_coeffs

    cp_upper_bound = []
    cp_lower_bound = []
    
    for temp in T:
        if temp <= trange[1]:
            cp_upper_bound.append(nasa7_cp(temp, upper_coeffs_low) * 8.314)
            cp_lower_bound.append(nasa7_cp(temp, lower_coeffs_low) * 8.314)
        else:
            cp_upper_bound.append(nasa7_cp(temp, upper_coeffs_high) * 8.314)
            cp_lower_bound.append(nasa7_cp(temp, lower_coeffs_high) * 8.314)
            
    return np.array(cp_upper_bound), np.array(cp_lower_bound)

def plot_cp_comparison(yaml_orig, yaml_opt, pkl_file_path, species_list):
    species_orig = load_species_data(yaml_orig)
    species_opt = load_species_data(yaml_opt)
    uncertainty_data = load_uncertainty_data(pkl_file_path)

    outdir = "opt_cp_curves_compare"
    os.makedirs(outdir, exist_ok=True)

    Tmin = 200
    Tmax = 3000
    T = np.linspace(Tmin, Tmax, 300)

    if uncertainty_data is None:
        print("Uncertainty data could not be loaded. Plotting without uncertainty bands.")

    for sp in species_list:
        cp_orig = get_cp_curve(species_orig, sp, T)
        cp_opt = get_cp_curve(species_opt, sp, T)

        if cp_orig is None:
            print(f"Species {sp} not found in original file.")
            continue
        if cp_opt is None:
            print(f"Species {sp} not found in optimized file.")
            continue

        plt.figure()
        plt.plot(T, cp_orig, 'k-', linewidth=2, label=f"{sp} (Original)")

        if uncertainty_data:
            cp_upper_bound, cp_lower_bound = get_uncertainty_band_curves(species_orig, uncertainty_data, sp, T)
            if cp_upper_bound is not None and cp_lower_bound is not None:
                plt.fill_between(T, cp_lower_bound, cp_upper_bound, color='gray', alpha=0.3, label="Uncertainty Band")

        plt.plot(T, cp_opt, 'r--', linewidth=2, label=f"{sp} (Optimized)")

        plt.xlabel("Temperature [K]")
        plt.ylabel("Cp [J/mol-K]")
        plt.title(f"Cp vs T for {sp}")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{sp}_Cp_compare_with_uncertainty.png"), dpi=300)
        plt.close()

if __name__ == "__main__":
    yaml_file_orig = "/data/TEST-THERMO-sens/FFCM1.0/FFCM1.yaml"
    yaml_file_opt = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/GA_ANALYSIS_with_opt_zeta/opt_itr_255885/optimized_results/OPTIMIZED_MECH/new_mech.yaml"
    pkl_file_path = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/unsrt.pkl"
    
    species_to_plot = ["AR", "H2", "O2", "H2O", "C", "CO", "OH*", "O", "HCO", "H"]
    plot_cp_comparison(yaml_file_orig, yaml_file_opt, pkl_file_path, species_to_plot)
    
'''
#####################################################


# usr this code to plot origional and optimized cp curves along with the uncertainty curves (MUQ-SAC) paper style direct run 

import numpy as np
import matplotlib.pyplot as plt
from ruamel.yaml import YAML
import os
import pickle

def nasa7_cp(T, coeffs):
    """
    Compute Cp/R from NASA7 coefficients at temperature T
    Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T**4
    """
    coeffs = np.array(coeffs)
    a1, a2, a3, a4, a5, a6, a7 = coeffs
    return a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4

def load_species_data(yaml_path):
    yaml = YAML(typ='safe')
    with open(yaml_path, 'r') as f:
        data = yaml.load(f)
    return data['species']

def get_cp_curve(species_data, sp, T):
    """Return Cp curve for a given species over T range"""
    for entry in species_data:
        if entry['name'].lower() == sp.lower():
            trange = entry['thermo']['temperature-ranges']
            coeffs_low = np.array(entry['thermo']['data'][0])
            coeffs_high = np.array(entry['thermo']['data'][1])

            cp_vals = []
            for temp in T:
                if temp <= trange[1]:
                    cp = nasa7_cp(temp, coeffs_low)
                else:
                    cp = nasa7_cp(temp, coeffs_high)
                cp_vals.append(cp * 8.314)
            return np.array(cp_vals)
    return None

def load_uncertainty_data(pkl_path):
    """Load uncertainty data from the pkl file"""
    if not os.path.exists(pkl_path):
        print(f"Error: The file '{pkl_path}' was not found.")
        return None
    try:
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
        return data
    except pickle.UnpicklingError as e:
        print(f"Error: Failed to unpickle the file. Details: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while loading the pkl file: {e}")
        return None
'''
def get_uncertainty_band_curves(species_orig_data, uncertainty_data, sp, T):
    """
    Calculates the upper and lower Cp curves based on uncertainty data.
    """
    orig_entry = None
    for entry in species_orig_data:
        if entry['name'].lower() == sp.lower():
            orig_entry = entry
            break
    
    if not orig_entry:
        print(f"Species {sp} not found in original data.")
        return None, None

    trange = orig_entry['thermo']['temperature-ranges']
    coeffs_low_orig = np.array(orig_entry['thermo']['data'][0])
    coeffs_high_orig = np.array(orig_entry['thermo']['data'][1])
     # 🌟 ADDED LINES TO PRINT NOMINAL COEFFICIENTS (p0) 🌟
    print(f"\n--- Nominal Coefficients for {sp} ---")
    print(f"Low-T Range (p0): {coeffs_low_orig}")
    print(f"High-T Range (p0): {coeffs_high_orig}")
    print("-------------------------------------------\n")

    try:
        nominal_obj_low = uncertainty_data[f"{sp}:Low"]
        nominal_obj_high = uncertainty_data[f"{sp}:High"]
        print(sp)
        print(sp)
    
        # 🌟 ADD THESE LINES TO PRINT THE VALUES INSTEAD OF THE OBJECT MEMORY LOCATION
        print(f"Low Object Attributes: ")
        print(f"  COV shape: {getattr(nominal_obj_low, 'cov', 'N/A').shape}")
        print(f"  Zeta Max x (Vector): {getattr(getattr(nominal_obj_low, 'zeta_max', None), 'x', 'N/A')}")
    
        print(f"High Object Attributes: ")
        print(f"  COV shape: {getattr(nominal_obj_high, 'cov', 'N/A').shape}")
        print(f"  Zeta Max x (Vector): {getattr(getattr(nominal_obj_high, 'zeta_max', None), 'x', 'N/A')}")
        # -------
        
        cov_low = getattr(nominal_obj_low, 'cov', None)
        zeta_max_obj_low = getattr(nominal_obj_low, 'zeta_max', None)
        zeta_max_vec_low = getattr(zeta_max_obj_low, 'x', None)
        
        cov_high = getattr(nominal_obj_high, 'cov', None)
        zeta_max_obj_high = getattr(nominal_obj_high, 'zeta_max', None)
        zeta_max_vec_high = getattr(zeta_max_obj_high, 'x', None)
    except (KeyError, AttributeError):
        print(f"Uncertainty data for {sp} is incomplete or not in expected format.")
        return None, None

    delta_low_coeffs = cov_low @ zeta_max_vec_low
    delta_high_coeffs = cov_high @ zeta_max_vec_high

    upper_coeffs_low = coeffs_low_orig.copy()
    lower_coeffs_low = coeffs_low_orig.copy()
    upper_coeffs_high = coeffs_high_orig.copy()
    lower_coeffs_high = coeffs_high_orig.copy()

    upper_coeffs_low[0:5] += delta_low_coeffs
    lower_coeffs_low[0:5] -= delta_low_coeffs
    upper_coeffs_high[0:5] += delta_high_coeffs
    lower_coeffs_high[0:5] -= delta_high_coeffs

    cp_upper_bound = []
    cp_lower_bound = []
    
    for temp in T:
        if temp <= trange[1]:
            cp_upper_bound.append(nasa7_cp(temp, upper_coeffs_low) * 8.314)
            cp_lower_bound.append(nasa7_cp(temp, lower_coeffs_low) * 8.314)
        else:
            cp_upper_bound.append(nasa7_cp(temp, upper_coeffs_high) * 8.314)
            cp_lower_bound.append(nasa7_cp(temp, lower_coeffs_high) * 8.314)
        
    return np.array(cp_upper_bound), np.array(cp_lower_bound)
'''



import numpy as np

def get_uncertainty_band_curves(species_orig_data, uncertainty_data, sp, T):
    """
    Calculates the upper and lower Cp curves based on uncertainty data, 
    and appends the nominal values, covariance, and zeta_max data for the species
    to a single master output file.
    """
    # -------------------------------------------------------------
    # 🌟 MODIFIED: Set a single output filename and use append mode 🌟
    # -------------------------------------------------------------
    MASTER_FILENAME = "all_species_uncertainty_data.txt"
    
    # 1. Search for original data
    orig_entry = None
    for entry in species_orig_data:
        if entry['name'].lower() == sp.lower():
            orig_entry = entry
            break
    
    if not orig_entry:
        print(f"Species {sp} not found in original data.")
        return None, None

    trange = orig_entry['thermo']['temperature-ranges']
    coeffs_low_orig = np.array(orig_entry['thermo']['data'][0])
    coeffs_high_orig = np.array(orig_entry['thermo']['data'][1])
    
    # 2. Extract Uncertainty Data
    try:
        nominal_obj_low = uncertainty_data[f"{sp}:Low"]
        nominal_obj_high = uncertainty_data[f"{sp}:High"]
        
        cov_low = getattr(nominal_obj_low, 'cov', 'N/A')
        zeta_max_obj_low = getattr(nominal_obj_low, 'zeta_max', None)
        zeta_max_vec_low = getattr(zeta_max_obj_low, 'x', 'N/A')
        
        cov_high = getattr(nominal_obj_high, 'cov', 'N/A')
        zeta_max_obj_high = getattr(nominal_obj_high, 'zeta_max', None)
        zeta_max_vec_high = getattr(zeta_max_obj_high, 'x', 'N/A')

    except (KeyError, AttributeError) as e:
        print(f"Uncertainty data for {sp} is incomplete or not in expected format. Error: {e}")
        return None, None
    
    # 3. Write data for THIS species to the master file in APPEND mode ('a')
    with open(MASTER_FILENAME, 'a') as f:
        f.write("\n" + "="*80 + "\n")
        f.write(f"UNCERTAINTY DATA FOR SPECIES: {sp}\n")
        f.write(f"Temperature Range: {trange[0]}K to {trange[2]}K (T_mid: {trange[1]}K)\n")
        f.write("="*80 + "\n\n")

        
        # --- LOW-TEMPERATURE RANGE DATA ---
        f.write(f"--- {sp} Low-Temperature Range ---\n")
        
        f.write(f"1. Nominal Coefficients (p0, Low-T):\n")
        np.savetxt(f, coeffs_low_orig.reshape(1, -1), fmt='%.12e', delimiter=', ')
        
        f.write(f"\n2. Covariance Matrix (Cov, Low-T):\n")
        if isinstance(cov_low, np.ndarray):
            np.savetxt(f, cov_low, fmt='%.12e', delimiter=', ')
        else:
            f.write(str(cov_low))
        
        f.write(f"\n3. Maximum Perturbation Vector (zeta_max, Low-T):\n")
        if isinstance(zeta_max_vec_low, np.ndarray):
            np.savetxt(f, zeta_max_vec_low.reshape(1, -1), fmt='%.12e', delimiter=', ')
        else:
            f.write(str(zeta_max_vec_low))
        f.write("\n")


        # --- HIGH-TEMPERATURE RANGE DATA ---
        f.write(f"--- {sp} High-Temperature Range ---\n")

        f.write(f"\n1. Nominal Coefficients (p0, High-T):\n")
        np.savetxt(f, coeffs_high_orig.reshape(1, -1), fmt='%.12e', delimiter=', ')

        f.write(f"\n2. Covariance Matrix (Cov, High-T):\n")
        if isinstance(cov_high, np.ndarray):
            np.savetxt(f, cov_high, fmt='%.12e', delimiter=', ')
        else:
            f.write(str(cov_high))

        f.write(f"\n3. Maximum Perturbation Vector (zeta_max, High-T):\n")
        if isinstance(zeta_max_vec_high, np.ndarray):
            np.savetxt(f, zeta_max_vec_high.reshape(1, -1), fmt='%.12e', delimiter=', ')
        else:
            f.write(str(zeta_max_vec_high))
        f.write("\n")
        
    print(f"✅ Appended uncertainty data for {sp} to: {MASTER_FILENAME}")
    # -------------------------------------------------------------
    # 🌟 END MODIFIED FILE SAVING 🌟
    # -------------------------------------------------------------

    # 4. Calculation Logic (Remains the same)
    delta_low_coeffs = cov_low @ zeta_max_vec_low
    delta_high_coeffs = cov_high @ zeta_max_vec_high

    upper_coeffs_low = coeffs_low_orig.copy()
    lower_coeffs_low = coeffs_low_orig.copy()
    upper_coeffs_high = coeffs_high_orig.copy()
    lower_coeffs_high = coeffs_high_orig.copy()

    upper_coeffs_low[0:5] += delta_low_coeffs
    lower_coeffs_low[0:5] -= delta_low_coeffs
    upper_coeffs_high[0:5] += delta_high_coeffs
    lower_coeffs_high[0:5] -= delta_high_coeffs

    cp_upper_bound = []
    cp_lower_bound = []
    
    # Assuming nasa7_cp and 8.314 are defined elsewhere and work as expected
    for temp in T:
        if temp <= trange[1]:
            cp_upper_bound.append(nasa7_cp(temp, upper_coeffs_low) * 8.314)
            cp_lower_bound.append(nasa7_cp(temp, lower_coeffs_low) * 8.314)
        else:
            cp_upper_bound.append(nasa7_cp(temp, upper_coeffs_high) * 8.314)
            cp_lower_bound.append(nasa7_cp(temp, lower_coeffs_high) * 8.314)
        
    return np.array(cp_upper_bound), np.array(cp_lower_bound)







def plot_cp_comparison(yaml_orig, yaml_opt, pkl_file_path, species_list):
    species_orig = load_species_data(yaml_orig)
    species_opt = load_species_data(yaml_opt)
    uncertainty_data = load_uncertainty_data(pkl_file_path)

    outdir = "opt_cp_curves_compare"
    os.makedirs(outdir, exist_ok=True)

    Tmin = 200
    Tmax = 3000
    T = np.linspace(Tmin, Tmax, 300)

    if uncertainty_data is None:
        print("Uncertainty data could not be loaded. Plotting without uncertainty bands.")

    for sp in species_list:
        cp_orig = get_cp_curve(species_orig, sp, T)
        cp_opt = get_cp_curve(species_opt, sp, T)

        if cp_orig is None:
            print(f"Species {sp} not found in original file.")
            continue
        if cp_opt is None:
            print(f"Species {sp} not found in optimized file.")
            continue

        plt.figure()
        
        # Plot in log scale, remove grids, and keep background empty
        plt.yscale('log')
        plt.grid(False)
        ax = plt.gca()
        ax.set_facecolor('white')

        plt.plot(T, cp_orig, 'k-', linewidth=2, label=f"{sp} (Original)")

        if uncertainty_data:
            cp_upper_bound, cp_lower_bound = get_uncertainty_band_curves(species_orig, uncertainty_data, sp, T)
            if cp_upper_bound is not None and cp_lower_bound is not None:
                # Plotting uncertainty limit curves without shade
                plt.plot(T, cp_upper_bound, 'g-', label="Uncertainty Limit (+)")
                plt.plot(T, cp_lower_bound, 'g-', label="Uncertainty Limit (-)")

        plt.plot(T, cp_opt, 'r-', linewidth=2, label=f"{sp} (Optimized)")

        plt.xlabel("Temperature [K]")
        plt.ylabel("Cp [J/mol-K]")
        plt.title(f"Cp vs T for {sp}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{sp}_Cp_compare_with_uncertainty.png"), dpi=300)
        plt.close()

if __name__ == "__main__":
    yaml_file_orig = "/data/TEST-THERMO-sens/FFCM1.0/FFCM1.yaml"
    yaml_file_opt = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/final_optimized_ga_400449/OPTIMIZED_MECH/new_mech.yaml"
    pkl_file_path = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/unsrt.pkl"
    
    species_to_plot = ["AR", "H2", "O2", "H2O", "C", "CO", "OH*", "O", "HCO", "H"]
    plot_cp_comparison(yaml_file_orig, yaml_file_opt, pkl_file_path, species_to_plot)





# usr this code to plot origional and optimized cp curves along with the uncertainty curves (MUQ-SAC) paper style module to be used in ga_opt_run
'''
import numpy as np
import matplotlib.pyplot as plt
from ruamel.yaml import YAML
import os
import pickle

def nasa7_cp(T, coeffs):
    """
    Compute Cp/R from NASA7 coefficients at temperature T
    Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T**4
    """
    coeffs = np.array(coeffs)
    a1, a2, a3, a4, a5, a6, a7 = coeffs
    return a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4

def load_species_data(yaml_path):
    """Load species thermo data from a YAML file."""
    yaml = YAML(typ='safe')
    if not os.path.exists(yaml_path):
        print(f"Error: The file '{yaml_path}' was not found.")
        return None
    with open(yaml_path, 'r') as f:
        data = yaml.load(f)
    return data['species']

def get_cp_curve(species_data, sp, T):
    """Return Cp curve for a given species over a temperature range."""
    if species_data is None:
        return None
        
    for entry in species_data:
        if entry['name'].lower() == sp.lower():
            trange = entry['thermo']['temperature-ranges']
            coeffs_low = np.array(entry['thermo']['data'][0])
            coeffs_high = np.array(entry['thermo']['data'][1])

            cp_vals = []
            for temp in T:
                if temp <= trange[1]:
                    cp = nasa7_cp(temp, coeffs_low)
                else:
                    cp = nasa7_cp(temp, coeffs_high)
                cp_vals.append(cp * 8.314)
            return np.array(cp_vals)
    return None

def load_uncertainty_data(pkl_path):
    """Load uncertainty data from the pkl file."""
    if not os.path.exists(pkl_path):
        print(f"Error: The file '{pkl_path}' was not found.")
        return None
    try:
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
        return data
    except pickle.UnpicklingError as e:
        print(f"Error: Failed to unpickle the file. Details: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while loading the pkl file: {e}")
        return None

def get_uncertainty_band_curves(species_orig_data, uncertainty_data, sp, T):
    """
    Calculates the upper and lower Cp curves based on uncertainty data.
    """
    orig_entry = None
    for entry in species_orig_data:
        if entry['name'].lower() == sp.lower():
            orig_entry = entry
            break
    
    if not orig_entry:
        print(f"Species {sp} not found in original data.")
        return None, None

    trange = orig_entry['thermo']['temperature-ranges']
    coeffs_low_orig = np.array(orig_entry['thermo']['data'][0])
    coeffs_high_orig = np.array(orig_entry['thermo']['data'][1])

    try:
        nominal_obj_low = uncertainty_data[f"{sp}:Low"]
        nominal_obj_high = uncertainty_data[f"{sp}:High"]
        
        zeta_max_vec_low = getattr(getattr(nominal_obj_low, 'zeta_max', None), 'x', None)
        zeta_max_vec_high = getattr(getattr(nominal_obj_high, 'zeta_max', None), 'x', None)
        
        cov_low = getattr(nominal_obj_low, 'cov', None)
        cov_high = getattr(nominal_obj_high, 'cov', None)
        
    except (KeyError, AttributeError):
        print(f"Uncertainty data for {sp} is incomplete or not in expected format.")
        return None, None
    
    if zeta_max_vec_low is None or zeta_max_vec_high is None or cov_low is None or cov_high is None:
        print(f"Uncertainty data for {sp} is missing key components.")
        return None, None

    delta_low_coeffs = cov_low @ zeta_max_vec_low
    delta_high_coeffs = cov_high @ zeta_max_vec_high

    upper_coeffs_low = coeffs_low_orig.copy()
    lower_coeffs_low = coeffs_low_orig.copy()
    upper_coeffs_high = coeffs_high_orig.copy()
    lower_coeffs_high = coeffs_high_orig.copy()

    upper_coeffs_low[0:5] += delta_low_coeffs
    lower_coeffs_low[0:5] -= delta_low_coeffs
    upper_coeffs_high[0:5] += delta_high_coeffs
    lower_coeffs_high[0:5] -= delta_high_coeffs

    cp_upper_bound = []
    cp_lower_bound = []
    
    for temp in T:
        if temp <= trange[1]:
            cp_upper_bound.append(nasa7_cp(temp, upper_coeffs_low) * 8.314)
            cp_lower_bound.append(nasa7_cp(temp, lower_coeffs_low) * 8.314)
        else:
            cp_upper_bound.append(nasa7_cp(temp, upper_coeffs_high) * 8.314)
            cp_lower_bound.append(nasa7_cp(temp, lower_coeffs_high) * 8.314)
        
    return np.array(cp_upper_bound), np.array(cp_lower_bound)

def plot_cp_comparison(yaml_orig, yaml_opt, pkl_file_path, species_list, output_dir):
    """
    Plots the Cp comparison for a list of species, including uncertainty bands.
    """
    species_orig = load_species_data(yaml_orig)
    species_opt = load_species_data(yaml_opt)
    uncertainty_data = load_uncertainty_data(pkl_file_path)

    os.makedirs(output_dir, exist_ok=True)

    Tmin = 200
    Tmax = 3000
    T = np.linspace(Tmin, Tmax, 300)

    if uncertainty_data is None:
        print("Uncertainty data could not be loaded. Plotting without uncertainty bands.")

    for sp in species_list:
        cp_orig = get_cp_curve(species_orig, sp, T)
        cp_opt = get_cp_curve(species_opt, sp, T)

        if cp_orig is None or cp_opt is None:
            print(f"Skipping species {sp} as data was not found in both files.")
            continue

        plt.figure()
        
        plt.yscale('log')
        plt.grid(False)
        ax = plt.gca()
        ax.set_facecolor('white')

        plt.plot(T, cp_orig, 'k-', linewidth=2, label=f"{sp} (Original)")

        if uncertainty_data:
            cp_upper_bound, cp_lower_bound = get_uncertainty_band_curves(species_orig, uncertainty_data, sp, T)
            if cp_upper_bound is not None and cp_lower_bound is not None:
                plt.plot(T, cp_upper_bound, 'g--', label="Uncertainty Limit (+)")
                plt.plot(T, cp_lower_bound, 'g--', label="Uncertainty Limit (-)")

        plt.plot(T, cp_opt, 'r--', linewidth=2, label=f"{sp} (Optimized)")

        plt.xlabel("Temperature [K]")
        plt.ylabel("Cp [J/mol-K]")
        plt.title(f"Cp vs T for {sp}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sp}_Cp_compare_with_uncertainty.png"), dpi=300)
        plt.close()
        print(f"✅ Plot for {sp} saved to {os.path.join(output_dir, f'{sp}_Cp_compare_with_uncertainty.png')}")


'''




