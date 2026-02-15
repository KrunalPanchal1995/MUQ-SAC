import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def match_curves(x1, y1, x2, y2):
    """
    Match two curves by interpolating the second curve onto the time points of the first.

    Args:
        x1, y1: Arrays representing the reference curve (time and OH profile).
        x2, y2: Arrays representing the simulated curve (time and OH profile).

    Returns:
        error_function: Mean squared error between the reference and simulated curves.
        y2_interpolated: Interpolated OH profile for the second curve.
    """
    # Interpolate y2 onto x1
    interpolator = interp1d(x2, y2, kind='linear', fill_value='extrapolate')
    y2_interpolated = interpolator(x1)
    
    # Calculate the error function (mean squared error)
    error_function = np.mean((y1 - y2_interpolated) ** 2)
    
    return error_function, y2_interpolated

def find_time(time, X):
    # Initialize dictionaries to store results
    times_below_X_max = {}
    X_below_max = {}
    times_after_X_max = {}
    X_after_max = {}
    list1 = [0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9]
    list2 = [0.99, 0.98, 0.97, 0.95, 0.91, 0.88]
           
    # Calculate max time
    max_X_idx = X.argmax()
    max_X_time = time[max_X_idx]
    X_at_max = X[max_X_idx]
        
   # Times and X_OH below dt-max at specific fractions
    for frac in list1:
        perc = int(frac * 100)
        X_target = frac * X_at_max
        idx = np.abs(X[:max_X_idx] - X_target).argmin()
        times_below_X_max[perc] = time[idx]
        X_below_max[perc] = X[idx]

    # Times and X_OH after max OH at specific fractions
    for frac in list2:
        perc = int(frac * 100)
        X_target = frac * X_at_max
        idx = np.abs(X[max_X_idx:] - X_target).argmin()
        times_after_X_max[perc] = time[max_X_idx+idx]
        X_after_max[perc] = X[max_X_idx+idx]
            
    #all_times = list(times_below_X_max.values()) + list(times_after_X_max.values())
    #all_X_OH = list(X_below_max.values()) + list(X_after_max.values())
    
     # Combine the results into lists
    time_results = (list(times_below_X_max.values())+ [max_X_time] + list(times_after_X_max.values()))
    X_results = (list(X_below_max.values())+ [X_at_max]+ list(X_after_max.values()))

    return time_results, X_results

    
# Example usage
if __name__ == "__main__":
    # Example synthetic data
    x1 = np.linspace(0, 10, 100)  # Reference time points
    y1 = np.exp(-0.5 * x1) * np.sin(2 * np.pi * x1)  # Reference OH profile
    
    x2 = np.linspace(0, 10, 80)  # Simulated time points
    y2 = 0.95 * np.exp(-0.5 * x2) * np.sin(2 * np.pi * x2)  # Simulated OH profile

    # Match curves and calculate error
    error, y2_matched = match_curves(x1, y1, x2, y2)

    # Print error function
    print(f"Mean Squared Error: {error}")

    # Plot the original and matched curves
    plt.figure(figsize=(10, 6))
    plt.plot(x1, y1, label="Reference OH Profile (x1, y1)", linewidth=2)
    plt.plot(x2, y2, label="Simulated OH Profile (x2, y2)", linestyle="--", alpha=0.7)
    plt.plot(x1, y2_matched, label="Interpolated Simulated Profile", linestyle="-.", alpha=0.9)
    plt.xlabel("Time")
    plt.ylabel("OH Profile")
    plt.title("Curve Matching - OH Profiles")
    plt.legend()
    plt.grid()
    plt.show()

