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
