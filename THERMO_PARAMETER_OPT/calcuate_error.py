

## this scripts takes two csv files as inputs( origional nominal sim , opt nominal sim) and compare the results (eror, absolute , relative , percentage)





'''
import os
import pandas as pd

def compute_mean_relative_error(original_csv, optimized_csv):
    try:
        df_orig = pd.read_csv(original_csv)
        df_red = pd.read_csv(reduced_csv)

        if len(df_orig) != len(df_red):
            print(f"⚠️ Row count mismatch: {original_csv}")
            return None

        # Match column names written earlier
        #orig_vals = df_orig["Obs(us)"].astype(float)
        orig_vals = df_orig["Obs(us)"].astype(float)
        red_vals = df_red["Reduced"].astype(float)

        rel_error = abs(red_vals - orig_vals) / orig_vals * 100
        return rel_error.mean()
    except Exception as e:
        print(f"❌ Error processing {original_csv}: {e}")
        return None

def evaluate_single_reduction_case(
    reduced_mech_folder,
    original_dataset_folder,
    target_types=("Tig", "RCM", "JSR", "Fls"),
    verbose=False
):
    """
    Compare reduced results to original and return mean relative error (%)
    Args:
        reduced_mech_folder (str): path to 'reduced_mechanism_xxx' folder
        original_dataset_folder (str): path to 'plot/Dataset'
        target_types (tuple): folders to look under (Tig, RCM, etc.)
        verbose (bool): whether to print per-case info
    Returns:
        float or None: mean relative error (%) for this reduction
    """
    
    reduced_plot_dir = os.path.join(reduced_mech_folder, "plot", "Dataset")
    if not os.path.isdir(reduced_plot_dir):
        print(f"❌ Reduced plot folder not found in {reduced_mech_folder}")
        return None

    all_errors = []

    for target_type in target_types:
        orig_dir = os.path.join(original_dataset_folder, target_type)
        red_dir = os.path.join(reduced_plot_dir, target_type)

        if not os.path.isdir(orig_dir) or not os.path.isdir(red_dir):
            continue

        for fname in os.listdir(orig_dir):
            if not fname.endswith(".csv"):#/home/user/Desktop/Mechanism_Reduction/trial_run_mech_reduction/test_opti_2/nsga2_results/reduced_0_8_83_457279/plot/Dataset/Tig/all_tig_data.csv
                continue
                                                                                        
            orig_path = os.path.join(orig_dir, fname)
            red_path = os.path.join(red_dir, fname)

            if not os.path.isfile(red_path):
                if verbose:
                    print(f"⚠️ Missing reduced file: {fname}")
                continue

            error = compute_mean_relative_error(orig_path, red_path)
            if error is not None:
                all_errors.append(error)
                if verbose:
                    print(f"{fname}: {error:.2f}%")

    if all_errors:
        mean_error = sum(all_errors) / len(all_errors)
        if verbose:
            print(f"\n✅ Mean Error for {os.path.basename(reduced_mech_folder)}: {mean_error:.2f}%")
        return mean_error
    else:
        print(f"⚠️ No valid CSV comparisons found in {reduced_mech_folder}")
        return None
        
        
'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def compare_and_plot_errors(original_csv, optimized_csv, output_plot_path="error_plot.png"):
    try:
        df_orig = pd.read_csv(original_csv)
        df_red = pd.read_csv(optimized_csv)

        if len(df_orig) != len(df_red):
            print(f"⚠️ Row count mismatch: {original_csv}")
            return None

        obs_vals = df_orig["Obs(us)"].astype(float)
        nominal_orig = df_orig["Nominal"].astype(float)
        nominal_red = df_red["Nominal"].astype(float)

        abs_error_orig = np.abs(obs_vals - nominal_orig)
        abs_error_red = np.abs(obs_vals - nominal_red)

        # Plot
        x = np.arange(1, len(obs_vals) + 1)  # 1, 2, 3, ..., N

        plt.figure(figsize=(10, 5))
        plt.scatter(x, abs_error_orig, label="Original Abs Error", color='blue', marker='o')
        plt.scatter(x, abs_error_red, label="optimized Abs Error", color='red', marker='x')
        plt.xlabel("Target Index")
        plt.ylabel("Absolute Error (us)")
        plt.title("Absolute Error Comparison")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()

        plt.savefig(output_plot_path)
        plt.close()
        print(f"✅ Plot saved to {output_plot_path}")

        # Compute statistics
        mae_orig = abs_error_orig.mean()
        mae_red = abs_error_red.mean()

        mre_orig = (abs_error_orig / obs_vals).mean() * 100
        mre_red = (abs_error_red / obs_vals).mean() * 100

        max_residual_orig = abs_error_orig.max()
        max_residual_red = abs_error_red.max()

        print("\n📊 Error Statistics:")
        print(f"🔵 Original: MAE = {mae_orig:.3f}, MRE = {mre_orig:.2f}%, Max Residual = {max_residual_orig:.2f}")
        print(f"🔴 optimized : MAE = {mae_red:.3f}, MRE = {mre_red:.2f}%, Max Residual = {max_residual_red:.2f}")

        return mae_orig, mae_red, mre_orig, mre_red, max_residual_orig, max_residual_red

    except Exception as e:
        print(f"❌ Error comparing and plotting: {e}")
        return None



if __name__ == "__main__":
    original_csv = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/GA_ANALYSIS_with_opt_zeta/origional_mech_results/plot/Dataset/Tig/all_tig_data.csv"
    optimized_csv = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/GA_ANALYSIS_with_opt_zeta/plot/Dataset/Tig/all_tig_data.csv"
    output_plot = "abs_error_plot.png"

    result = compare_and_plot_errors(original_csv, optimized_csv, output_plot)

    if result:
        mae_o, mae_r, mre_o, mre_r, max_o, max_r = result

        print("\n📈 Summary of Error Metrics:")
        print(f"🔹 Mean Absolute Error (Original): {mae_o:.3f} us")
        print(f"🔸 Mean Absolute Error (Reduced) : {mae_r:.3f} us")

        print(f"🔹 Mean Relative Error (Original): {mre_o:.2f}%")
        print(f"🔸 Mean Relative Error (Reduced) : {mre_r:.2f}%")

        print(f"🔹 Max Residual Error (Original) : {max_o:.2f} us")
        print(f"🔸 Max Residual Error (Reduced)  : {max_r:.2f} us")
    else:
        print("❌ Failed to compute or plot errors.")


























