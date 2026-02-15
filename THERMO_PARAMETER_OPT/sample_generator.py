import os
import numpy as np
import matplotlib.pyplot as plt

def sample_plot_all_species(unsrt_data, design_matrix, outdir="sample_plots", seed=None):
    """
    Same high‑level purpose as before, but now the species order is
    taken from the insertion order of `unsrt_data`, guaranteeing that
    blocks of 10 ζ‑values align with (species, Low/High) pairs.

    Every other step (θ‑matrices, envelope bands, plotting, etc.) is
    unchanged.
    """
    # ---------------------------------------------------------------------
    # 0. bookkeeping & input checks
    # ---------------------------------------------------------------------
    os.makedirs(outdir, exist_ok=True)

    # ── preserve order ───────────────────────────────────────────────────
    species_order = []
    _seen = set()
    for key in unsrt_data:                # dicts are ordered in Py ≥3.7
        sp = key.split(":")[0]
        if sp not in _seen:
            species_order.append(sp)
            _seen.add(sp)

    n_sp = len(species_order)
    expected_cols = n_sp * 10             # 5 ζ_low + 5 ζ_high per species
    '''
    if design_matrix.shape[1] != expected_cols:
        raise ValueError(
            f"design_matrix has {design_matrix.shape[1]} columns but "
            f"{n_sp} species ⇒ need {expected_cols}."
        ) '''

    rng = np.random.default_rng(seed)
    #if design_matrix.shape[0] < 100:
        #raise ValueError("Need at least 100 rows in design_matrix.")
    #chosen_rows = rng.choice(design_matrix.shape[0], size=100, replace=False)
    chosen_rows = rng.choice(len(design_matrix), size=100, replace=False)

    # ---------------------------------------------------------------------
    # 1. pre‑compute θ‑matrices and cache stuff that never changes
    # ---------------------------------------------------------------------
    R   = 8.314
    T_lo = np.linspace(300, 1000, 100)
    T_hi = np.linspace(1000, 3000, 100)

    θ_lo = np.vstack([T_lo/T_lo, T_lo, T_lo**2, T_lo**3, T_lo**4]).T
    θ_hi = np.vstack([T_hi/T_hi, T_hi, T_hi**2, T_hi**3, T_hi**4]).T

    cache = {}
    for sp in species_order:
        lo = unsrt_data[f"{sp}:Low"]
        hi = unsrt_data[f"{sp}:High"]

        n_lo = np.asarray(lo.NominalParams, dtype=float).flatten()
        n_hi = np.asarray(hi.NominalParams, dtype=float).flatten()
        print("species name : \t ", sp)
        print("low nominal \t ", n_lo)
        print("high nominal ", n_hi)
        

        # “Envelope” (max/min) curves
        cp_lo_max = n_lo + lo.cov @ lo.zeta_max.x
        cp_lo_min = n_lo - lo.cov @ lo.zeta_max.x
        cp_hi_max = n_hi + hi.cov @ hi.zeta_max.x
        cp_hi_min = n_hi - hi.cov @ hi.zeta_max.x

        cache[sp] = dict(
            nom_lo    = R * (θ_lo @ n_lo),
            nom_hi    = R * (θ_hi @ n_hi),
            env_lo_hi = R * (θ_lo @ cp_lo_max),
            env_lo_lo = R * (θ_lo @ cp_lo_min),
            env_hi_hi = R * (θ_hi @ cp_hi_max),
            env_hi_lo = R * (θ_hi @ cp_hi_min),
            cov_lo    = lo.cov,
            cov_hi    = hi.cov,
        )

    # ---------------------------------------------------------------------
    # 2. loop over 100 random rows × all species
    # ---------------------------------------------------------------------
    for row_idx in chosen_rows:
        row = design_matrix[row_idx]

        for sp_idx, sp in enumerate(species_order):
            base = sp_idx * 10
            ζ_lo = row[base      : base + 5].reshape(-1, 1)
            ζ_hi = row[base + 5  : base + 10].reshape(-1, 1)

            c = cache[sp]

            cp_lo = R * (θ_lo @ (c["cov_lo"] @ ζ_lo)).flatten() + c["nom_lo"]
            cp_hi = R * (θ_hi @ (c["cov_hi"] @ ζ_hi)).flatten() + c["nom_hi"]

            # ── plotting ────────────────────────────────────────────────
            fig, ax = plt.subplots(figsize=(8, 5))

            ax.plot(T_lo, cp_lo, label="Sample (<1000 K)",  color="red")
            ax.plot(T_hi, cp_hi, label="Sample (>1000 K)", color="green")

            ax.plot(T_lo, c["nom_lo"], color="blue", label="Nominal")
            ax.plot(T_hi, c["nom_hi"], color="blue")

            ax.plot(T_lo, c["env_lo_hi"], "k--", label="Uncert. limits")
            ax.plot(T_lo, c["env_lo_lo"], "k--")
            ax.plot(T_hi, c["env_hi_hi"], "k--")
            ax.plot(T_hi, c["env_hi_lo"], "k--")

            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel(r"$C_p$ (J mol$^{-1}$ K$^{-1}$)")
            ax.set_ylim(0.92 * c["env_lo_lo"].min(),
                        1.15 * c["env_hi_hi"].max())
            ax.legend(fontsize=9)
            ax.grid(False)

            plt.savefig(f"{outdir}/{sp}_{row_idx}.pdf", bbox_inches="tight")
            plt.close(fig)

    print(f"✓ 100 rows × {n_sp} species ⇒ {len(chosen_rows)*n_sp} Cp‑plots "
          f"in '{os.path.abspath(outdir)}'")



















def generate_samples_csv(unsrt_data, design_matrix, csv_path="samples.csv"):
    """
    For *every* row of `design_matrix`, convert the ζ‑values into
    sampled NASA‑polynomial coefficients

        a_sample = a_nominal + cov @ ζ

    and write the results to `csv_path`.

    ── Column order ──────────────────────────────────────────
    For each species in the insertion order of `unsrt_data`
    (same as in your plotting script) the 10 columns are:

        {sp}:Low:a1  …  {sp}:Low:a5   {sp}:High:a1  …  {sp}:High:a5

    So the CSV has `n_rows` rows × (`n_species` × 10) columns.
    """
    # ── 1. work out the canonical species ordering ───────────
    species_order, _seen = [], set()
    for key in unsrt_data:                         # dicts are ordered (Py ≥ 3.7)
        sp = key.split(":")[0]
        if sp not in _seen:
            species_order.append(sp)
            _seen.add(sp)

    n_sp = len(species_order)
    expected_cols = n_sp * 10
    if design_matrix.shape[1] != expected_cols:
        raise ValueError(
            f"design_matrix has {design_matrix.shape[1]} columns but "
            f"expect {expected_cols} (= {n_sp} species × 10)."
        )

    # ── 2. build a column header list ────────────────────────
    header = []
    for sp in species_order:
        header += [f"{sp}:Low:a{i}"  for i in range(1, 6)]
        header += [f"{sp}:High:a{i}" for i in range(1, 6)]

    # ── 3. pre‑pull the pieces that never change ─────────────
    cache = {}
    for sp in species_order:
        lo = unsrt_data[f"{sp}:Low"]
        hi = unsrt_data[f"{sp}:High"]
        cache[sp] = dict(
            n_lo=np.asarray(lo.NominalParams, dtype=float).flatten(),   # (5,)
            n_hi=np.asarray(hi.NominalParams, dtype=float).flatten(),
            cov_lo=np.asarray(lo.cov, dtype=float),                    # (5,5)
            cov_hi=np.asarray(hi.cov, dtype=float),
        )

    # ── 4. iterate over every row of the design matrix ───────
    rows_out = []
    for row_idx in range(design_matrix.shape[0]):
        row = design_matrix[row_idx]
        sample_row = []

        for sp_idx, sp in enumerate(species_order):
            base = sp_idx * 10
            ζ_lo = row[base     : base + 5].reshape(-1, 1)  # (5,1)
            ζ_hi = row[base + 5 : base + 10].reshape(-1, 1) # (5,1)

            c = cache[sp]

            # a_sample = a_nom + cov @ ζ
            a_lo = c["n_lo"] + (c["cov_lo"] @ ζ_lo).flatten()
            a_hi = c["n_hi"] + (c["cov_hi"] @ ζ_hi).flatten()

            sample_row.extend(a_lo.tolist() + a_hi.tolist())

        rows_out.append(sample_row)

    # ── 5. write the CSV ─────────────────────────────────────
    df = pd.DataFrame(rows_out, columns=header)
    df.to_csv(csv_path, index=False)
    print(f"✅ Wrote {df.shape[0]} rows × {df.shape[1]} cols to '{os.path.abspath(csv_path)}'")















