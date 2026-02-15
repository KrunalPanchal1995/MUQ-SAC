# use the code to manually check the cp values for any row of DM actual parameters
import numpy as np
import matplotlib.pyplot as plt

# Universal Gas Constant (R) in J/(mol*K)
R_universal = 8.314

# --- 1. Paste Your 10 Coefficient Values Here ---
# First 5 are for the LOW temperature range: a1_low, a2_low, a3_low, a4_low, a5_low
# Next 5 are for the HIGH temperature range: a1_high, a2_high, a3_high, a4_high, a5_high

# Placeholder coefficients - REPLACE THESE WITH YOUR ACTUAL 10 VALUES
# Example placeholder values for CO2 (approximate for demonstration)
# Coeffs a1 to a5 for T < 1000 K (Low range)


#,,3.0286393350130427,0.00809669021629345,-1.9654955110854693e-05,2.024579690304718e-08,-7.351936560333525e-12,3.0031319026152032,0.0019479075334090695,-7.220450249613295e-07,1.4659353961180686e-10,-1.1353167509449764e-14,2.7208123353384455,-0.003493340329644302,6.841746693828702e-06,-6.133657260966047e-09,2.1081697290664115e-12,2.8107746252993535,-0.0013881657925994737,7.900258355780593e-07,-1.8413039304880014e-10,1.5226892393845127e-14,3.5550028849976343,-0.0031616494279557953,9.947005177196775e-06,-9.704728698731253e-09,3.2364392601321546e-12,3.393439592195521,0.0005935670319979884,-1.361661412417531e-07,2.2924878925176343e-11,-1.6961662373770917e-15,3.309608421500987,-0.0023556920215019126,6.6202016387717935e-06,-5.570832038714221e-09,1.749089689171592e-12,1.8698273669389058,0.002524360735452354,-7.444312651623818e-07,1.0893972412688786e-10,-6.320872127503075e-15,3.414084336037936,-0.0006720550359469652,1.069408293321507e-06,8.804371674031281e-10,-9.118733341640031e-13,3.0342083654731993,0.001011356537334395,-3.0046225712345586e-07,3.6094177838295234e-11,-1.1953968708977788e-15,2.581278995521631,-0.0002933643520526034,7.112330084206271e-07,-7.317421880535051e-10,2.6687922810790787e-13,2.639656828990138,-0.00019640369191418526,1.0636993410171547e-07,-1.6118044663431955e-11,7.796654298157511e-16,4.046730471836387,-0.0035190153738692007,1.403223427530646e-05,-1.3418558673823958e-08,4.375073841217457e-12,3.6686903962699295,0.0023982743665579647,-6.460559519757179e-07,1.028369648405518e-10,-7.28123502564134e-15,3.145499076157819,0.0004684181215027582,-1.9559500975714462e-06,3.149842986132058e-09,-1.3648475300018464e-12,2.4590785276799236,0.0013433911424961336,-4.227421674758367e-07,6.726221530878297e-11,-4.027161789788238e-15,2.956603749584837,0.0004912630286778383,-3.9205621826526406e-07,8.402294129020453e-12,6.120858694083e-15,3.1755960329886466,-0.00019566233698692145,1.1584801566128484e-07,-2.7764897689881817e-11,2.3168988473198857e-15


coeffs_Low = np.array([2.893360598630726,0.00042321929949450256,-3.3775339964672984e-07,7.238511403964778e-12,5.273072422702961e-15]) 
# Coeffs a1 to a5 for T > 1000 K (High range)
coeffs_High = np.array([3.082020756970549,-0.00016856158994892805,9.980216945686259e-08,-2.391924461010923e-11,1.9959904367348083e-15]) 
# ----------------------------------------------------

# --- 2. Define Temperature Ranges ---
T_low_start, T_low_end = 300, 1000  # Low temperature range [300 1000] K
T_high_start, T_high_end = 1000, 1800  # High temperature range [1000 1800] K

# Create a smooth array of temperature points for plotting
T_low = np.linspace(T_low_start, T_low_end, 50)
T_high = np.linspace(T_high_start, T_high_end, 50)


# --- 3. Function to Calculate Cp/R ---
def calculate_cp_over_R(T, coeffs):
    """Calculates Cp/R using the NASA 5-term polynomial (a1 to a5)."""
    a1, a2, a3, a4, a5 = coeffs
    
    # Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    cp_over_R = (a1 + 
                 a2 * T + 
                 a3 * T**2 + 
                 a4 * T**3 + 
                 a5 * T**4)
    return cp_over_R

# --- 4. Calculate Cp Values ---

# Calculate Cp/R for the low-temperature range
cp_over_R_low = calculate_cp_over_R(T_low, coeffs_low)
Cp_low = cp_over_R_low * R_universal

# Calculate Cp/R for the high-temperature range
cp_over_R_high = calculate_cp_over_R(T_high, coeffs_high)
Cp_high = cp_over_R_high * R_universal

# --- 5. Plotting ---

plt.figure(figsize=(10, 6))

# Plot Low Temperature Curve
plt.plot(T_low, Cp_low, label=f'Low Temp. Range ({T_low_start}-{T_low_end} K)', 
         color='blue', linestyle='-', linewidth=2)

# Plot High Temperature Curve
plt.plot(T_high, Cp_high, label=f'High Temp. Range ({T_high_start}-{T_high_end} K)', 
         color='red', linestyle='--', linewidth=2)

# Add a marker at the crossover point (1000 K)
# This uses the calculated value at the end of the low range and start of the high range
Cp_at_1000_low = Cp_low[-1]
Cp_at_1000_high = Cp_high[0]
plt.plot(T_low_end, Cp_at_1000_low, 'o', color='blue', markersize=7, label=f'$C_p$ at {T_low_end} K (Low)')
plt.plot(T_high_start, Cp_at_1000_high, 'x', color='red', markersize=7, label=f'$C_p$ at {T_high_start} K (High)')

plt.title('Specific Heat Capacity ($C_p$) vs. Temperature ($T$)')
plt.xlabel('Temperature (K)')
plt.ylabel(f'Specific Heat Capacity, $C_p$ (J/mol·K) [using R={R_universal} J/mol·K]')
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.show()

print(f"\nCalculated $C_p$ at {T_low_end} K (from Low Coeffs): {Cp_at_1000_low:.2f} J/mol·K")
print(f"Calculated $C_p$ at {T_high_start} K (from High Coeffs): {Cp_at_1000_high:.2f} J/mol·K")
if abs(Cp_at_1000_low - Cp_at_1000_high) == 0.0:
    print("\nNote: The $C_p$ values at the transition temperature (1000 K) are NOT continuous with the placeholder data.")
else:
    print("\nNote: The $C_p$ values are continuous at the transition temperature (1000 K).")
