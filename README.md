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

```bash
conda env create -f env/environment.yml
conda activate rbfe
---

## Running hydration free energy
---
```bash
python rbfe/run_hydration_absolute.py

---

Outputs per-seed results and MBAR QC plots under:
---
```bash
results/hydration_abs/
results/qc/

## Running complex absolute free energy
---
```bash
python rbfe/run_complex_absolute.py

