# openmm-rbfe
Open-source RBFE (free energy perturbation FEP) using OpenMM, perses, and OpenFF for small molecule drug discovery

# OpenMM RBFE Side Project

This repository contains a project focused on building an open-source relative binding free energy (RBFE) workflow for small-molecule
drug discovery using OpenMM, perses, and OpenFF.

## Project Goal
To benchmark RBFE (FEP) predictions on a congeneric ligand series against public experimental binding affinity data and to establish quality-control criteria for identifying reliable and unreliable free energy estimates.

## Methods
- Molecular dynamics: OpenMM
- Alchemical free energy: perses
- Small-molecule force field: OpenFF
- Analysis: MBAR (pymbar), MDTraj, RDKit

## Planned Deliverables
- Reproducible RBFE workflow (protein + ligand transformations)
- ΔΔG ranking and comparison with experimental data
- Analysis of convergence, uncertainty, and failure cases
- Decision-oriented summary for medicinal chemistry prioritization

