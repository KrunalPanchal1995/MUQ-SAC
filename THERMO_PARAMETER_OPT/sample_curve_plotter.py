'''
import os
import numpy as np
import matplotlib.pyplot as plt

def sample_plot_all_species(unsrt_data, design_matrix, outdir="sample_plots"):
    """
    Iterate over every row in `design_matrix`. Each row has 100 entries arranged
    as 10 blocks of 10 (5 low-T zeta, 5 high-T zeta) – one block per species.
    For every (row, species) pair, generate a Cp–T curve PDF that shows
    ▸ the nominal Cp,
    ▸ the sample Cp curve, and
    ▸ ±uncertainty‐envelope curves.

    Parameters
    ----------
    unsrt_data : dict-like
        Keys look like "{species}:Low" and "{species}:High"; each value has
        attributes `.cov`, `.NominalParams`, and `.zeta_max.x`.
    design_matrix : (n_samples, 100) ndarray-like
        Latin-hyper-cube (or similar) samples.  For a given row, columns
        [0:10) are species-0 (5 low, 5 high), columns [10:20) species-1, etc.
    outdir : str, optional
        Destination directory for PDFs.

    Notes
    -----
    • The species order is inferred from `unsrt_data` and assumed to match the
      block order in each row of `design_matrix`.
    • Uses a fixed low-T grid (300–1000 K) and high-T grid (1000–3000 K), each
      with 100 points.
    """

    # ---------- 0. housekeeping ------------------------------------------------
    os.makedirs(outdir, exist_ok=True)

    # ---------- 1. discover species & cache constant pieces --------------------
    species_names = sorted({key.split(":")[0] for key in unsrt_data})

    T_low   = np.linspace(300, 1000, 100)
    T_high  = np.linspace(1000, 3000, 100)
    Tscale  = 1.0

    theta_low  = np.vstack([ T_low/T_low,
                             T_low/Tscale,
                             (T_low/Tscale)**2,
                             (T_low/Tscale)**3,
                             (T_low/Tscale)**4 ]).T        # (100,5)
    theta_high = np.vstack([ T_high/T_high,
                             T_high/Tscale,
                             (T_high/Tscale)**2,
                             (T_high/Tscale)**3,
                             (T_high/Tscale)**4 ]).T        # (100,5)

    R = 8.314

    # Cache nominal & ±uncertainty envelopes for every species
    cache = {}
    for sp in species_names:
        lo  = unsrt_data[f"{sp}:Low"]
        hi  = unsrt_data[f"{sp}:High"]

        npar_lo  = np.asarray(lo.NominalParams).flatten()
        npar_hi  = np.asarray(hi.NominalParams).flatten()

        # envelopes (vector length 5 each)
        cp_max_lo  = npar_lo + lo.cov @ lo.zeta_max.x
        cp_min_lo  = npar_lo - lo.cov @ lo.zeta_max.x
        cp_max_hi  = npar_hi + hi.cov @ hi.zeta_max.x
        cp_min_hi  = npar_hi - hi.cov @ hi.zeta_max.x

        # store everything that never changes across rows
        cache[sp] = dict(
            theta_lo   = theta_low,
            theta_hi   = theta_high,
            nom_lo     = R * theta_low  @ npar_lo,
            nom_hi     = R * theta_high @ npar_hi,
            env_lo_max = R * theta_low  @ cp_max_lo,
            env_lo_min = R * theta_low  @ cp_min_lo,
            env_hi_max = R * theta_high @ cp_max_hi,
            env_hi_min = R * theta_high @ cp_min_hi,
        )

    # ---------- 2. iterate rows & species --------------------------------------
    for row_idx, row in enumerate(design_matrix):
        for sp_idx, sp in enumerate(species_names):
            # slice the 10-col block for this species
            base = sp_idx * 10
            zeta_lo  = np.asarray(row[base : base + 5]).reshape(-1, 1)  # (5,1)
            zeta_hi  = np.asarray(row[base + 5 : base + 10]).reshape(-1, 1)

            c          = cache[sp]         # shortcut
            samp_lo_cp = R * (c["theta_lo"]  @ (unsrt_data[f"{sp}:Low"].cov  @ zeta_lo)).flatten() + c["nom_lo"]
            samp_hi_cp = R * (c["theta_hi"]  @ (unsrt_data[f"{sp}:High"].cov @ zeta_hi)).flatten() + c["nom_hi"]

            # ---------------------- 3. plotting --------------------------------
            fig, ax = plt.subplots(figsize=(10, 6))

            ax.plot(T_low,  samp_lo_cp,  label="Sample (T < 1000 K)", color="red",   lw=0.8)
            ax.plot(T_high, samp_hi_cp,  label="Sample (T > 1000 K)", color="green", lw=0.8)

            ax.plot(T_low,  c["nom_lo"], label="Nominal",            color="blue")
            ax.plot(T_high, c["nom_hi"], color="blue")

            ax.plot(T_low,  c["env_lo_max"], "k--", label="Unsrt. limits")
            ax.plot(T_low,  c["env_lo_min"], "k--")
            ax.plot(T_high, c["env_hi_max"], "k--")
            ax.plot(T_high, c["env_hi_min"], "k--")

            ax.set_xlabel("Temperature (K)", fontsize=20)
            ax.set_ylabel(r"Heat Capacity $C_p$ (J mol$^{-1}$ K$^{-1}$)", fontsize=20)
            ax.set_ylim(0.92 * c["env_lo_min"].min(),
                        1.15 * c["env_hi_max"].max())
            ax.legend(fontsize=12, loc="best")

            for spine in ax.spines.values():
                spine.set_edgecolor("black")
                spine.set_linewidth(1)

            ax.grid(False)

            fname = f"{outdir}/{sp}_row{row_idx:04d}.pdf"
            plt.savefig(fname, bbox_inches="tight")
            plt.close(fig)

    print(f"All Cp curves have been written to → {os.path.abspath(outdir)}")

'''

import os
import numpy as np
import matplotlib.pyplot as plt

def sample_plot_all_species(unsrt_data, design_matrix, outdir="sample_plots", seed=None):
    """
    From `design_matrix` choose 100 random rows, then for *each* of the 10 species
    build and save a Cp(T) curve PDF.  Output files are named
        <species>_<rowindex>.pdf
    where <rowindex> is the original row number in `design_matrix`.

    Parameters
    ----------
    unsrt_data : dict-like
        Must contain keys like "H2:Low", "H2:High", each with .cov,
        .NominalParams, and .zeta_max.x.
    design_matrix : (N, 100) ndarray
        Latin-hypercube (or similar) samples: 10 blocks × 10 ζ-values each.
    outdir : str, optional
        Destination directory for PDFs (created if absent).
    seed : int or None, optional
        Give a seed for reproducible random-row selection.
    """
    # -------------------------------------------------------------------------
    # 0. set-up & random row selection
    # -------------------------------------------------------------------------
    os.makedirs(outdir, exist_ok=True)

    rng = np.random.default_rng(seed)
    n_rows = design_matrix.shape[0]
    if n_rows < 100:
        raise ValueError(f"design_matrix has only {n_rows} rows; need ≥100.")
    chosen_rows = rng.choice(n_rows, size=100, replace=False)


    # -------------------------------------------------------------------------
    # 1. discover the 10 species (order inferred from unsrt_data keys)
    # -------------------------------------------------------------------------
    species = {k.split(":")[0] for k in unsrt_data}  # len == 10
    #species = [k for k in unsrt_data]
    n_sp = len(species)
    print("Sample curve plotter species lenght", n_sp)
    #if n_sp != 10:
        #raise RuntimeError(f"Expected 10 species but found {n_sp}: {species}")

    # -------------------------------------------------------------------------
    # 2. pre-compute θ-matrices and cache everything that never changes across
    #    rows – avoids doing this inside the nested loops
    # -------------------------------------------------------------------------
    R = 8.314
    T_lo  = np.linspace(300, 1000, 100)
    T_hi  = np.linspace(1000, 3000, 100)

    theta_lo = np.vstack([T_lo/T_lo, T_lo, T_lo**2, T_lo**3, T_lo**4]).T
    theta_hi = np.vstack([T_hi/T_hi, T_hi, T_hi**2, T_hi**3, T_hi**4]).T

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
        )

    # -------------------------------------------------------------------------
    # 3. main loop: random rows × 10 species
    # -------------------------------------------------------------------------
    for row_idx in chosen_rows:
        row = design_matrix[row_idx]

        for sp_idx, sp in enumerate(species):
            base = sp_idx * 10
            z_lo = row[base : base+5].reshape(-1, 1)
            z_hi = row[base+5 : base+10].reshape(-1, 1)
            #z_hi = row[base : base+5].reshape(-1, 1)
            #z_lo = row[base+5 : base+10].reshape(-1, 1)

            c = cache[sp]

            sample_lo = R * (theta_lo @ (c["cov_lo"] @ z_lo)).flatten() + c["nom_lo"]
            sample_hi = R * (theta_hi @ (c["cov_hi"] @ z_hi)).flatten() + c["nom_hi"]

            # ---------- plotting ------------------------------------------------
            fig, ax = plt.subplots(figsize=(8, 5))

            ax.plot(T_lo,  sample_lo, label="Sample (<1000 K)",  color="red")
            ax.plot(T_hi,  sample_hi, label="Sample (>1000 K)", color="green")

            ax.plot(T_lo,  c["nom_lo"], color="blue", label="Nominal")
            ax.plot(T_hi,  c["nom_hi"], color="blue")

            ax.plot(T_lo,  c["env_lo_hi"], "k--", label="Uncert. limits")
            ax.plot(T_lo,  c["env_lo_lo"], "k--")
            ax.plot(T_hi,  c["env_hi_hi"], "k--")
            ax.plot(T_hi,  c["env_hi_lo"], "k--")

            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel(r"$C_p$ (J mol$^{-1}$ K$^{-1}$)")
            ax.set_ylim(0.92 * c["env_lo_lo"].min(),
                        1.15 * c["env_hi_hi"].max())
            ax.legend(fontsize=9)
            ax.grid(False)

            fname = f"{outdir}/{sp}_{row_idx}.pdf"
            plt.savefig(fname, bbox_inches="tight")
            plt.close(fig)

    print(f"✓ 100 rows × 10 species ⇒ {len(chosen_rows)*10} Cp-plots in '{os.path.abspath(outdir)}'")



