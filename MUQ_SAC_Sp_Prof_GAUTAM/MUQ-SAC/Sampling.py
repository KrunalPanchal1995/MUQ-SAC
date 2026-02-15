import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
list1 = [0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9]
list2 = [0.99, 0.98, 0.97, 0.95, 0.91, 0.88]

class Sampling():
    def __init__(self, time, X):
        self.time = time
        self.X = X
        self.dt_max_times = []
        self.X_at_dt_max_list = []
        self.onset_times = []
        self.X_at_onset_times = []
        self.max_X_times = []
        self.X_at_max_list = []
        self.times_below_X_max_list = {int(frac * 100): [] for frac in list1}
        self.X_before_X_max_list = {int(frac * 100): [] for frac in list1}
        self.times_after_X_max_list = {int(frac * 100): [] for frac in list2}
        self.X_after_X_max_list = {int(frac * 100): [] for frac in list2}
        
    def calculate_max_X_time(self):
        max_X_idx = self.X.argmax()
        max_X_time = self.time[max_X_idx]
        X_at_max = self.X[max_X_idx]
        self.max_X_times.append(max_X_time)
        self.X_at_max_list.append(X_at_max)
        return max_X_time, X_at_max

    def calculate_below_max(self, max_X_time, X_at_max):
        max_X_idx = self.X.argmax()
        max_X_time = self.time[max_X_idx]
        X_at_max = self.X[max_X_idx]
        for frac in list1:
            perc = int(frac * 100)
            X_target = frac * X_at_max
            idx = np.abs(self.X[:max_X_idx] - X_target).argmin()
            self.times_below_X_max_list[perc] = self.time[idx]
            self.X_before_X_max_list[perc] = self.X[idx]
        return self.times_below_X_max_list, self.X_before_X_max_list

    def calculate_after_max(self, max_X_time, X_at_max):
        max_X_idx = self.X.argmax()
        max_X_time = self.time[max_X_idx]
        X_at_max = self.X[max_X_idx]
        for frac in list2:
            perc = int(frac * 100)
            X_target = frac * X_at_max
            idx = np.abs(self.X[max_X_idx:] - X_target).argmin()
            self.times_after_X_max_list[perc]=self.time[max_X_idx+idx]
            self.X_after_X_max_list[perc] = self.X[max_X_idx+idx]
        return self.times_after_X_max_list, self.X_after_X_max_list
        
    def generate_csv_files(
        times_below_max, X_below_max, max_X_times,
        X_at_max_list, times_after_max_list, X_after_max_list,
        key):
        # Create dictionaries
        time_dict = {}
        time_dict.update({
            f'times_below_max_{perc}': times_below_max[perc]
            for perc in times_below_max
        })

        concentration_dict = {}
        concentration_dict.update({
            f'X_below_max_{perc}': X_below_max[perc]
            for perc in X_below_max
        })
        time_dict.update({
            'time_max': max_X_times
        })
        
        concentration_dict.update({
            'max_X': X_at_max_list
        })

        # Add times after max to the dictionaries
        time_dict.update({
            f'times_after_max_{perc}': times_after_max_list[perc]
            for perc in times_after_max_list
        })

        concentration_dict.update({
            f'X_after_max_{perc}': X_after_max_list[perc]
            for perc in X_after_max_list
        })
    
        # Convert to DataFrames
        time_df = pd.DataFrame(time_dict)
        concentration_df = pd.DataFrame(concentration_dict)
        
        # Generate correct column names dynamically
        desired_time_columns = [f'times_below_max_{int(frac * 100)}' for frac in list1] + ['time_max'] + [f'times_after_max_{int(frac * 100)}' for frac in list2]
        time_df = time_df[desired_time_columns]
        # Save to CSV files
        time_file_path = f"time_values_{key}.csv"
        concentration_file_path = f"concentration_values_{key}.csv"
        time_df.to_csv(time_file_path, index=False)
        concentration_df.to_csv(concentration_file_path, index=False)
