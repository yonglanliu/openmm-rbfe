from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

QC_RE = re.compile(r"^(?P<ligand>.+)_seed(?P<seed>\d+)_mbar_qc\.csv$")

def find_ddg_summary(results_dir: Path) -> Path | None:
    p = results_dir / "ddg_summary_repeats.csv"
    return p if p.exists() else None

def load_ddg_summary(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)

def scan_qc_csv(results_dir: Path) -> list[dict]:
    qc_dir = results_dir / "qc"
    if not qc_dir.exists():
        return []
    items = []
    for p in sorted(qc_dir.glob("*_mbar_qc.csv")):
        m = QC_RE.match(p.name)
        if not m:
            continue
        items.append(
            {
                "ligand": m.group("ligand"),
                "seed": int(m.group("seed")),
                "path": p,
            }
        )
    return items

def load_qc_csv(path: Path) -> pd.DataFrame:
    # expected columns from save_mbar_qc():
    # i, le_i, ls_i, le_ip1, ls_ip1, dG_adj_kJmol, dG_adj_err_kJmol, dG_cum_to_i_kJmol, dG_cum_err_to_i_kJmol
    return pd.read_csv(path)

def available_ligands(qc_items: list[dict]) -> list[str]:
    return sorted({it["ligand"] for it in qc_items})

def available_seeds(qc_items: list[dict], ligand: str) -> list[int]:
    return sorted({it["seed"] for it in qc_items if it["ligand"] == ligand})

def qc_path_for(qc_items: list[dict], ligand: str, seed: int) -> Path | None:
    for it in qc_items:
        if it["ligand"] == ligand and it["seed"] == seed:
            return it["path"]
    return None

def summarize_ligand(lig_tag, qc_dir="results/qc", outdir="results"):
    """
    Run QC summarization programmatically and return (rows, mean, std, sem, out_csv)
    instead of printing to console.
    """
    import glob
    import numpy as np
    import os
    import csv

    pattern = os.path.join(qc_dir, f"{lig_tag}_complex_seed*_mbar_qc.csv")
    qc_files = sorted(glob.glob(pattern))

    if not qc_files:
        raise FileNotFoundError(f"No QC files found matching pattern: {pattern}")

    rows = []  # (seed, dG, sem)
    for path in qc_files:
        basename = os.path.basename(path)
        try:
            after_seed = basename.split("seed", 1)[1]
            seed_str = after_seed.split("_", 1)[0]
            seed = int(seed_str)
        except Exception:
            raise RuntimeError(f"Could not parse seed from QC filename: {basename}")

        # Reuse the existing function
        dG_cum, dG_cum_err = extract_final_d_G_from_qc(path)
        rows.append((seed, dG_cum, dG_cum_err))

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

    return rows, mean, std, sem, out_csv
