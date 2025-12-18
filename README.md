# OpenMM-based Absolute Free Energy Pipeline

This repository implements a **reproducible absolute free energy (AFE) pipeline**
based on **OpenMM**, **openmmtools**, and **MBAR**, supporting:

- Absolute hydration free energies (ΔG_hyd)
- Absolute protein–ligand complex decoupling
- Multi-seed uncertainty estimation
- MBAR-based QC and λ-space diagnostics
- Visualization dashboard (Streamlit)

---

## Why this repository exists

Absolute hydration free energy is a **canonical benchmark** for validating:
- force fields (GAFF / AMBER),
- alchemical protocols,
- λ schedules,
- sampling stability and uncertainty.

This pipeline is designed as a **calibration + validation step**
before running **absolute binding free energies**.

---

## Repository structure

openmm-rbfe/
├── rbfe/ # core simulation scripts
│ ├── run_hydration_absolute.py
│ ├── run_complex_absolute.py
│
├── tools/ # protein & ligand preprocessing
├── viz/ # Streamlit visualization app
├── env/ # conda environment
├── data/ # input data (examples only)
├── results/ # (ignored) simulation outputs


---

## Environment setup

### Local (macOS / CPU)
```bash
# from repo root
conda env create -f env/environment.yml
conda activate rbfe

python -c "import openmm; print('OpenMM OK')"
```
---
### HPC (GPU/CUDA)
```bash
# Use conda environment
git clone https://github.com/yonglanliu/openmm-rbfe.git
cd openmm-rbfe

# Create env
conda env create -f env/environment.yml
conda activate rbfe

# verify OpenMM platforms (should include CUDA on GPU node)
python - << 'PY'
import openmm as mm
print([mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])
PY
```

## Run: Absolute hydration free energy (ΔG_hyd)
Prepare input ligands:
Put ligand SDF files under data/ligands/
Update the ligand paths in rbfe/run_hydration_absolute.py (or use CLI if enabled)
Run:
```bash
python rbfe/run_hydration_absolute.py
```
Ouputs (typical):
* results/hydration_abs/<ligand>/seedXXXX/hydration_result.json
* results/qc/<ligand>_seedXXXX_mbar_qc.csv
* results/qc/<ligand>_seedXXXX_mbar_qc.png (+ cumulative plot)
---

## Run: Absolute complex decoupling (ΔG_complex)
Inputs:
* Protein PDB (processed; missing atoms fixed; keep desired ligand via resname)
* Ligand SDF file that matches the ligand in the PDB (same chemistry)
Run:
```bash
python rbfe/run_complex_absolute.py
```
Outputs (typical):
* results/complex_abs/<target>/<ligand>/seedXXXX/complex_result.json
* QC CSV/plots under results/qc/

### Minimal QC: What to check
Recommended quick checks:
* MBAR adjacent ΔG contributions: no single window dominates unexpectedly
* Cumulative ΔG curve: should be reasonably smooth
* Multi-seed repeatability: seed-to-seed std should be small relative to target precision