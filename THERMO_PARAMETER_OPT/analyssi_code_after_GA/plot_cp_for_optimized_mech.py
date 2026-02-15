## this script takes input the new yaml file after optimization and plots the cp curves for the species given by the user
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
'''
def plot_cp_for_species(yaml_path, species_list):
    species_data = load_species_data(yaml_path)

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

                plt.plot(T, cp_vals, label=sp)
        if not found:
            print(f"Species {sp} not found in file.")

    plt.xlabel("Temperature [K]")
    plt.ylabel("Cp [J/mol-K]")
    plt.title("Cp vs T for selected species")
    plt.legend()
    plt.grid(True)
    plt.show()
'''




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
    yaml_file = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/GA_ANALYSIS_with_opt_zeta/new_mech.yaml"
    species_to_plot = ["AR", "H2", "O2", "H2O" , "C" , "CO", "OH*", "O" , "HCO" , "H"]  # <-- manually change this list
    plot_cp_for_species(yaml_file, species_to_plot)

