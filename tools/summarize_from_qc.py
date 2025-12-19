"""Author: Yonglan Liu
Created: 2025-12
Project: OpenMM-based Free Energy Pipeline (Absolute Hydration)
"""

"""
Usage:
    python tools/summarize_complex_from_qc.py --ligand_tag MBN
"""

#!/usr/bin/env python

import os
import glob
import csv
import argparse
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize complex-leg free energies across seeds "
            "by reading MBAR QC CSV files."
        )
    )
    parser.add_argument(
        "--ligand_tag",
        "--lig_tag",
        dest="lig_tag",
        required=True,
        help=(
            "Ligand tag used in QC filenames, e.g. 'MBN' if files "
            "are like results/qc/MBN_complex_seed2025_mbar_qc.csv."
        ),
    )
    parser.add_argument(
        "--qc_dir",
        dest="qc_dir",
        default="results/qc",
        help="Directory containing *_mbar_qc.csv files (default: results/qc).",
    )
    parser.add_argument(
        "--outdir",
        dest="outdir",
        default="results",
        help="Directory to write summary CSV (default: results).",
    )
    return parser.parse_args()


import csv

def extract_final_dG_from_qc(csv_path):
    last_row = None
    with open(csv_path, "r") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        for row in reader:
            if not row:
                continue
            last_row = row

    if last_row is None:
        raise RuntimeError(f"No data rows found in QC file: {csv_path}")

    # Last two columns are cumulative dG and its error
    dG_cum = float(last_row[-2])
    dG_cum_err = float(last_row[-1])
    return dG_cum, dG_cum_err

def summarize_ligand(lig_tag, qc_dir="results/qc", outdir="results", leg="complex"):
    """
    Summarize QC for a given ligand and leg.
    leg: "complex" or "hydration" (used in QC filename pattern)
    Expects files like:
      results/qc/<LIG>_complex_seedXXXX_mbar_qc.csv
      results/qc/<LIG>_hydration_seedXXXX_mbar_qc.csv
    """
    import glob
    import numpy as np
    import os
    import csv

    # Note: leg is inserted into the filename
    pattern = os.path.join(qc_dir, f"{lig_tag}_{leg}_seed*_mbar_qc.csv")
    qc_files = sorted(glob.glob(pattern))

    if not qc_files:
        raise FileNotFoundError(f"No QC files found matching pattern: {pattern}")

    rows = []  # (seed, dG, sem)
    for path in qc_files:
        basename = os.path.basename(path)
        # parse seed from e.g. MBN_complex_seed2025_mbar_qc.csv
        try:
            after_seed = basename.split("seed", 1)[1]
            seed_str = after_seed.split("_", 1)[0]
            seed = int(seed_str)
        except Exception:
            raise RuntimeError(f"Could not parse seed from QC filename: {basename}")

        dG_cum, dG_cum_err = extract_final_dG_from_qc(path)
        rows.append((seed, dG_cum, dG_cum_err))

    dGs = np.array([r[1] for r in rows], dtype=float)
    mean = float(dGs.mean())
    std = float(dGs.std(ddof=1)) if len(dGs) > 1 else 0.0
    sem = float(std / np.sqrt(len(dGs))) if len(dGs) > 1 else 0.0

    os.makedirs(outdir, exist_ok=True)
    # Include leg in the output filename so complex vs hydration are separate
    out_csv = os.path.join(outdir, f"{lig_tag}_{leg}_summary_repeats_from_qc.csv")

    with open(out_csv, "w") as f:
        f.write("seed,dG_kJmol,sem_kJmol\n")
        for seed, dG, s in rows:
            f.write(f"{seed},{dG:.6f},{s:.6f}\n")
        f.write(f"MEAN,{mean:.6f},\n")
        f.write(f"STD,{std:.6f},\n")
        f.write(f"SEM,{sem:.6f},\n")

    return rows, mean, std, sem, out_csv

def main():
    args = parse_args()
    lig_tag = args.lig_tag
    qc_dir = args.qc_dir
    outdir = args.outdir

    pattern = os.path.join(qc_dir, f"{lig_tag}_complex_seed*_mbar_qc.csv")
    qc_files = sorted(glob.glob(pattern))

    if not qc_files:
        raise SystemExit(f"No QC files found matching pattern: {pattern}")

    print(f"Found {len(qc_files)} QC files for ligand tag '{lig_tag}':")
    for path in qc_files:
        print("  ", path)

    rows = []  # (seed, dG, sem)
    for path in qc_files:
        # Extract seed from filename, assuming ..._seed<SEED>_mbar_qc.csv
        basename = os.path.basename(path)
        # e.g. MBN_complex_seed2025_mbar_qc.csv
        # split on "seed" then on "_"
        try:
            after_seed = basename.split("seed", 1)[1]
            seed_str = after_seed.split("_", 1)[0]
            seed = int(seed_str)
        except Exception:
            raise RuntimeError(
                f"Could not parse seed from QC filename: {basename}"
            )

        dG_cum, dG_cum_err = extract_final_dG_from_qc(path)
        print(
            f"  seed {seed}: dG = {dG_cum:.3f} ± {dG_cum_err:.3f} kJ/mol "
            f"(from {basename})"
        )
        rows.append((seed, dG_cum, dG_cum_err))

    # Aggregate statistics across seeds (use dG values only)
    dGs = np.array([r[1] for r in rows], dtype=float)
    mean = float(dGs.mean())
    std = float(dGs.std(ddof=1)) if len(dGs) > 1 else 0.0
    sem = float(std / np.sqrt(len(dGs))) if len(dGs) > 1 else 0.0

    os.makedirs(outdir, exist_ok=True)
    out_csv = os.path.join(outdir, f"{lig_tag}_complex_summary_repeats_from_qc.csv")

    with open(out_csv, "w") as f:
        f.write("seed,dG_complex_kJmol,sem_kJmol\n")
        for seed, dG, s in rows:
            f.write(f"{seed},{dG:.6f},{s:.6f}\n")
        f.write(f"MEAN,{mean:.6f},\n")
        f.write(f"STD,{std:.6f},\n")
        f.write(f"SEM,{sem:.6f},\n")

    print("\n===== SUMMARY ACROSS SEEDS =====")
    for seed, dG, s in rows:
        print(f"seed {seed}: dG = {dG:.3f} ± {s:.3f} kJ/mol")
    print("-------------------------------")
    print(f"Mean dG = {mean:.3f} kJ/mol")
    print(f"Std dG  = {std:.3f} kJ/mol")
    print(f"SEM dG  = {sem:.3f} kJ/mol")
    print(f"\nWrote summary CSV: {out_csv}")


if __name__ == "__main__":
    main()


