import matplotlib.pyplot as plt
from MechanismParser import Parser
import numpy as np
import os,sys

Nominal_Mech = sys.argv[1]
Perutrbed_Mech = sys.argv[2]

nominal = Parser(Nominal_Mech).mech
perturbed = Parser(Perutrbed_Mech).mech

species = f"{sys.argv[3]}"
#print(nominal)
for index,spec_data in enumerate(nominal["species"]):
    if spec_data["name"] == species:
       nom_NASA_7 =  spec_data["thermo"]
for index,spec_data in enumerate(perturbed["species"]):
    if spec_data["name"] == species:
       per_NSAS_7 = spec_data["thermo"]

def cp_component_of_enthalpy(T,coeff):
    R = 8.314
    
    a1, a2, a3, a4, a5, a6, a7 = coeff
    
    H_RT = (
    a1 + a2 * T / 2 + a3 * T**2 / 3 +
    a4 * T**3 / 4 + a5 * T**4 / 5 
    )
    H_RT_1000 = H_RT   # Returns enthalpy H in cal/mol
    return H_RT_1000    

def get_Enthalpy_Curves(T,coeff):
    R = 8.314
    
    a1, a2, a3, a4, a5, a6, a7 = coeff
    
    H_RT = (
    a1 + a2 * T / 2 + a3 * T**2 / 3 +
    a4 * T**3 / 4 + a5 * T**4 / 5 + a6 / T
    )
    H_RT_1000 = H_RT*T   # Returns enthalpy H in cal/mol
    return H_RT_1000

def get_Thermo_Curves(T1,coeff):
    R = 8.314
    T = np.linspace(T1[0],T1[-1],50)
    
    a1, a2, a3, a4, a5, a6, a7 = coeff
    # Calculate enthalpy (H) at temperature T
    Cp_R = a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4
    
    Cp = Cp_R * R
    
    H_RT = (
    a1 + a2 * T / 2 + a3 * T**2 / 3 +
    a4 * T**3 / 4 + a5 * T**4 / 5 + a6 / T
    )
    H = H_RT * R * T  # Returns enthalpy H in cal/mol

    # Calculate entropy (S) at temperature T
    S_R = (
    a1 * np.log(T) + a2 * T + a3 * T**2 / 2 +
    a4 * T**3 / 3 + a5 * T**4 / 4 + a7
    )
    S = S_R * R  # Returns entropy S in
    return Cp,H, S

T = np.linspace(298,nom_NASA_7["temperature-ranges"][1],50)
T_ = np.linspace(nom_NASA_7["temperature-ranges"][1],nom_NASA_7["temperature-ranges"][2],50)


cp, h , s = get_Thermo_Curves(T,nom_NASA_7["data"][0])


cp_, h_ , s_ = get_Thermo_Curves(T,per_NSAS_7["data"][0])


cph, hh , sh = get_Thermo_Curves(T_,nom_NASA_7["data"][1])


cph_, hh_ , sh_ = get_Thermo_Curves(T_,per_NSAS_7["data"][1])


per_NSAS_7["data"][1][5] = get_Enthalpy_Curves(1000,per_NSAS_7["data"][0]) - 1000*cp_component_of_enthalpy(1000,per_NSAS_7["data"][1])

cph__, hh__ , sh__ = get_Thermo_Curves(T_,per_NSAS_7["data"][1])

fig = plt.figure()
plt.plot(T,cp,"k-", label="Nominal_low")
plt.plot(T_,cph,"k-.",label="Nominal_high")
plt.plot(T,cp_,"b-", label=f"Cp_low")
plt.plot(T_,cph_,"g-",label=f"Cp_high")
plt.legend()
plt.ylabel(rf"Heat Capacity (Cp / J/mol.K, {sys.argv[3]})")
plt.xlabel("Temperature (K)")
plt.savefig("Comparison_cp.pdf",bbox_inches="tight")

fig = plt.figure()
plt.plot(T,h,"k-", label="Nominal_low")
plt.plot(1000,8.314*get_Enthalpy_Curves(1000,per_NSAS_7["data"][0]),'o',label='H_o(1000K)')
plt.plot(T_,hh,"k-.",label="Nominal_high")

plt.plot(T,h_,"b-", label=f"H_low")
plt.plot(T_,hh_,"g-",label=f"H_high")
plt.plot(T_,hh__,"r--")
plt.legend()
plt.ylabel(rf"Enthalpy (J/mol.K, {sys.argv[3]})")
plt.xlabel("Temperature (K)")
plt.savefig("Comparison_h.pdf",bbox_inches="tight")

fig = plt.figure()
plt.plot(T,s,"k-", label="Nominal_low")
plt.plot(T_,sh,"k-.",label="Nominal_high")
plt.plot(T,s_,"b-", label=f"S_low")
plt.plot(T_,sh_,"g-",label=f"S_high")
plt.legend()
plt.ylabel(rf"Entropy (J/K, {sys.argv[3]})")
plt.xlabel("Temperature (K)")
plt.savefig("Comparison_s.pdf",bbox_inches="tight")

