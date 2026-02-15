
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# helper: load files
# -----------------------------------------------------------------------------
def load_unsrt_data(pkl_path):
    """Load the uncertainty data dictionary from pickle file."""
    with open(pkl_path, "rb") as f:
        return pickle.load(f)

def load_design_matrix(txt_path):
    """
    Load design matrix from text file where each row ends with a comma.
    Example row: 0.12,0.34,0.56,
    """
    rows = []
    with open(txt_path, "r") as f:
        for line in f:
            vals = line.strip().split(",")
            vals = [float(v) for v in vals if v != ""]
            if vals:
                rows.append(vals)
    return np.array(rows)

# -----------------------------------------------------------------------------
# helper: plotting
# -----------------------------------------------------------------------------
def plot_cp_curves(sp, row_idx, T_lo, T_hi, sample_lo, sample_hi, cache, outdir):
    """Plot Cp curves for one species and one row index."""
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(T_lo, sample_lo, label="Sample (<1000 K)",  color="red")
    ax.plot(T_hi, sample_hi, label="Sample (>1000 K)",  color="green")

    ax.plot(T_lo, cache["nom_lo"], color="blue", label="Nominal")
    ax.plot(T_hi, cache["nom_hi"], color="blue")

    ax.plot(T_lo, cache["env_lo_hi"], "k--", label="Uncert. limits")
    ax.plot(T_lo, cache["env_lo_lo"], "k--")
    ax.plot(T_hi, cache["env_hi_hi"], "k--")
    ax.plot(T_hi, cache["env_hi_lo"], "k--")

    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$C_p$ (J mol$^{-1}$ K$^{-1}$)")
    ax.set_ylim(0.92 * cache["env_lo_lo"].min(),
                1.15 * cache["env_hi_hi"].max())
    ax.legend(fontsize=9)
    ax.grid(False)

    fname = f"{outdir}/{sp}_{row_idx}.pdf"
    plt.savefig(fname, bbox_inches="tight")
    plt.close(fig)

# -----------------------------------------------------------------------------
# main logic
# -----------------------------------------------------------------------------
def sample_plot_all_species(pkl_path, design_matrix_path,
                            outdir="sample_plots",
                            seed=None,
                            param_outfile="opt_param.txt"):
    """
    From `design_matrix` choose 100 random rows, then for each species
    build and save Cp(T) curves. Also save the Cp parameter values into `opt_param.txt`.
    """
    # 0. set-up & random row selection
    os.makedirs(outdir, exist_ok=True)

    unsrt_data = load_unsrt_data(pkl_path)
    design_matrix = load_design_matrix(design_matrix_path)

    rng = np.random.default_rng(seed)
    n_rows = design_matrix.shape[0]
    chosen_rows = rng.choice(n_rows, size=1, replace=False)

    # 1. discover species
    #species = {k.split(":")[0] for k in unsrt_data}
    #print("Sample curve plotter species length:", len(species))
    species = sorted({k.split(":")[0] for k in unsrt_data})
    for sp in species:
        assert f"{sp}:Low" in unsrt_data and f"{sp}:High" in unsrt_data, f"Missing data for {sp}"
    print("Sample curve plotter species length:", len(species))

    # --- NEW: dump unsrt data (nominal params + cov) ---
    unsrt_dump_file = os.path.join(outdir, "unsrt_data_check.txt")
    with open(unsrt_dump_file, "w") as f_dump:
        for sp in species:
            for tag in ["Low", "High"]:
                key = f"{sp}:{tag}"
                entry = unsrt_data[key]

                f_dump.write(f"{sp}_{tag.lower()} - Nominal Params:\n")
                f_dump.write(",".join(map(str, np.asarray(entry.NominalParams).flatten())) + "\n")

                f_dump.write(f"{sp}_{tag.lower()} - Covariance Matrix:\n")
                np.savetxt(f_dump, entry.cov, fmt="%.6e", delimiter=",")
                f_dump.write("\n")  # blank line for readability
    print(f"✓ Unsrt data written to '{unsrt_dump_file}'")

    # 2. pre-compute θ-matrices
    R = 8.314
    T_lo  = np.linspace(300, 1000, 100)
    T_hi  = np.linspace(1000, 3000, 100)

    theta_lo = np.vstack([T_lo/T_lo, T_lo, T_lo**2, T_lo**3, T_lo**4]).T
    theta_hi = np.vstack([T_hi/T_hi, T_hi, T_hi**2, T_hi**3, T_hi**4]).T

    # Precompute cache
    cache = {}
    for sp in species:
        lo = unsrt_data[f"{sp}:Low"]
        hi = unsrt_data[f"{sp}:High"]

        n_lo = np.asarray(lo.NominalParams).flatten()
        n_hi = np.asarray(hi.NominalParams).flatten()

        cp_lo_max = n_lo + lo.cov @ lo.zeta_max.x
        cp_lo_min = n_lo - lo.cov @ lo.zeta_max.x
        cp_hi_max = n_hi + hi.cov @ hi.zeta_max.x
        cp_hi_min = n_hi - hi.cov @ hi.zeta_max.x

        cache[sp] = dict(
            nom_lo    = R * theta_lo @ n_lo,
            nom_hi    = R * theta_hi @ n_hi,
            env_lo_hi = R * theta_lo @ cp_lo_max,
            env_lo_lo = R * theta_lo @ cp_lo_min,
            env_hi_hi = R * theta_hi @ cp_hi_max,
            env_hi_lo = R * theta_hi @ cp_hi_min,
            cov_lo    = lo.cov,
            cov_hi    = hi.cov,
            n_lo      = n_lo,
            n_hi      = n_hi
        )

    # 3. open param_outfile
    with open(param_outfile, "w") as f_out:
        # main loop
        for row_idx in chosen_rows:
            row = design_matrix[row_idx]

            for sp_idx, sp in enumerate(species):
                base = sp_idx * 10
                z_lo = row[base : base+5].reshape(-1, 1)
                z_hi = row[base+5 : base+10].reshape(-1, 1)

                c = cache[sp]

                # Cp(T) samples
                sample_lo = R * (theta_lo @ (c["cov_lo"] @ z_lo)).flatten() + c["nom_lo"]
                sample_hi = R * (theta_hi @ (c["cov_hi"] @ z_hi)).flatten() + c["nom_hi"]

                # --- write parameters (Nominal + Cov*zeta) ---
                cp_params_lo = (c["n_lo"] + (c["cov_lo"] @ z_lo).flatten()).tolist()
                cp_params_hi = (c["n_hi"] + (c["cov_hi"] @ z_hi).flatten()).tolist()

                #all_params = cp_params_lo + cp_params_hi
                #f_out.write(f"{sp},{row_idx}," + ",".join(map(str, all_params)) + "\n")
                # --- write parameters in AR:Low_a1=... format ---
                for i, val in enumerate(cp_params_lo, start=1):
                    f_out.write(f"{sp}:Low_a{i}=\t{val}\n")
                for i, val in enumerate(cp_params_hi, start=1):
                    f_out.write(f"{sp}:High_a{i}=\t{val}\n")
                f_out.write("\n")  # blank line between species


                # --- plotting ---
                plot_cp_curves(sp, row_idx, T_lo, T_hi, sample_lo, sample_hi, c, outdir)

    print(f"✓ {len(chosen_rows)} rows × {len(species)} species ⇒ "
          f"{len(chosen_rows)*len(species)} Cp-plots in '{os.path.abspath(outdir)}'")
    print(f"✓ Cp parameter values written to '{param_outfile}'")

# -----------------------------------------------------------------------------
# entry point (no CLI args, just edit paths here)
# -----------------------------------------------------------------------------
def main():
    # ---- Put your paths here ----
    pkl_path = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/unsrt.pkl"
    design_matrix_path = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/GA_ANALYSIS_with_opt_params/solution_zeta_values.save"
    outdir = "opt_plots"        # or "/custom/path/sample_plots"
    param_outfile = os.path.join(outdir, "opt_param.txt")

    # Call main processing function
    sample_plot_all_species(
        pkl_path,
        design_matrix_path,
        outdir=outdir,
        seed=42,
        param_outfile=param_outfile
    )

if __name__ == "__main__":
    main()


