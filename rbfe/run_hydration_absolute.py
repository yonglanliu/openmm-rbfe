"""
Author: Yonglan Liu
Created: 2025-12
Project: OpenMM-based Free Energy Pipeline (Absolute Hydration)

What this script does
---------------------
This module computes *absolute hydration free energies* (ΔG_hyd) for small
molecules using OpenMM + openmmtools alchemy + MBAR (PyMBAR).

It is designed as a reproducible, pipeline-friendly entry point that:
  1) builds a solvated ligand system (GAFF + AMBER TIP3P),
  2) constructs an alchemical system via AbsoluteAlchemicalFactory,
  3) runs a two-stage decoupling schedule (electrostatics -> sterics),
  4) evaluates reduced potentials u_kln for all states (k,l,n),
  5) estimates ΔG and uncertainty using MBAR, and
  6) writes per-seed results and QC artifacts (CSV + plots).

Why hydration free energy?
--------------------------
Hydration free energy is a standard benchmark for validating:
  - force field choices (e.g., GAFF/SMIRNOFF),
  - alchemical protocol correctness (lambda naming, PME treatment),
  - numerical stability (NaNs, bad contacts),
  - sampling quality and overlap (MBAR diagnostics, seed sensitivity).

This is a recommended calibration step before running protein–ligand binding
free energies (absolute or relative).

Design principles
-----------------
- Reproducibility-first:
  Key parameters (schedule, temperature, friction, timestep, platform, seed,
  sampling lengths) are configurable and recorded to disk.
- Clear separation of concerns:
  system build -> relaxation -> alchemical sampling -> estimator -> QC output
- QC transparency:
  Window-wise ΔG contributions and cumulative curves are exported for quick
  diagnostics.

Typical outputs
---------------
- results/hydration_abs/<ligand>/seed<seed>/hydration_result.json
- results/qc/<ligand>_seed<seed>_mbar_qc.csv
- results/qc/<ligand>_seed<seed>_mbar_qc.png
- results/qc/<ligand>_seed<seed>_mbar_qc_cum.png
- results/ddg_summary_repeats.csv (if main() runs multiple seeds)

Run
---
    python run_hydration_absolute.py
"""

import os
import json
import numpy as np

from openmm import unit
import openmm as mm
import openmm.app as app

from rdkit import Chem
from rdkit.Chem import AllChem

from openff.toolkit import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator

from openmmtools.alchemy import AlchemicalRegion
from pymbar import MBAR


# -----------------------------
#          Ligand prep
# -----------------------------
def _rdkit_embed_sdf(in_sdf: str, out_sdf: str, random_seed: int = 2025) -> None:
    """Ensure ligand has 3D coords; write a new SDF with explicit H + 3D."""
    suppl = Chem.SDMolSupplier(in_sdf, removeHs=False)
    m = suppl[0]
    if m is None:
        raise ValueError(f"Failed to read SDF: {in_sdf}")

    m = Chem.AddHs(m, addCoords=True)
    if m.GetNumConformers() == 0:
        AllChem.EmbedMolecule(m, randomSeed=int(random_seed))
    AllChem.UFFOptimizeMolecule(m)

    w = Chem.SDWriter(out_sdf)
    w.write(m)
    w.close()


# ---------------------------------
# System build (GAFF + AMBER TIP3P)
# ---------------------------------
def build_solvated_system(ligand_sdf: str, padding_nm: float = 1.2, tmp_sdf: str | None = None):
    """
    Build solvated OpenMM System+Topology+Positions for a single ligand
    using GAFF (via AmberTools) + AMBER TIP3P water.

    Returns: (system, topology, positions)
    """
    if tmp_sdf is None:
        tmp_sdf = ligand_sdf.replace(".sdf", ".3d.sdf")
    _rdkit_embed_sdf(ligand_sdf, tmp_sdf)

    offmol = Molecule.from_file(tmp_sdf)
    lig_top = offmol.to_topology().to_openmm()
    positions = offmol.conformers[0].to_openmm()

    forcefield = app.ForceField("amber14/tip3p.xml")
    gaff = GAFFTemplateGenerator(molecules=[offmol], forcefield="gaff-2.11")
    forcefield.registerTemplateGenerator(gaff.generator)

    modeller = app.Modeller(lig_top, positions)
    modeller.addSolvent(forcefield, model="tip3p", padding=padding_nm * unit.nanometer)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.AllBonds,
        rigidWater=True,
    )
    return system, modeller.topology, modeller.positions


def get_ligand_atom_indices(topology: app.Topology) -> list[int]:
    """Heuristic: ligand is the only non-water residue."""
    water_names = {"HOH", "WAT", "SOL", "TIP3"}
    atom_indices = []
    for atom in topology.atoms():
        if atom.residue.name.upper() not in water_names:
            atom_indices.append(atom.index)
    if len(atom_indices) == 0:
        raise RuntimeError("Could not identify ligand atoms in topology.")
    return atom_indices


# -------------------------------------------------
# Alchemical system (openmmtools 0.25.3 compatible)
# -------------------------------------------------
def alchemical_system(reference_system: mm.System, alchemical_atoms: list[int]) -> mm.System:
    """Build alchemical system using openmmtools AbsoluteAlchemicalFactory (0.25.3 style)."""
    from openmmtools import alchemy

    region = AlchemicalRegion(alchemical_atoms=alchemical_atoms)

    factory = alchemy.AbsoluteAlchemicalFactory(
        consistent_exceptions=False,
        alchemical_pme_treatment="exact",
    )

    alch_system = factory.create_alchemical_system(
        reference_system=reference_system,
        alchemical_regions=region,
        alchemical_regions_interactions=frozenset(),
    )
    return alch_system


# -----------------------------
# Utilities: lambda parameters
# -----------------------------
def _lambda_param_names(system: mm.System) -> list[str]:
    """Return names of global parameters that look like lambda parameters."""
    names = set()
    for fi in range(system.getNumForces()):
        f = system.getForce(fi)
        if hasattr(f, "getNumGlobalParameters") and hasattr(f, "getGlobalParameterName"):
            for pi in range(f.getNumGlobalParameters()):
                n = f.getGlobalParameterName(pi)
                if "lambda" in n.lower():
                    names.add(n)
    return sorted(names)


def make_decouple_schedule(n_elec: int = 6, n_vdw: int = 12):
    """
    Two-stage decoupling schedule:
      Stage 1: lambda_electrostatics 1 -> 0 with lambda_sterics = 1
      Stage 2: lambda_sterics       1 -> 0 with lambda_electrostatics = 0

    Returns list of (lambda_electrostatics, lambda_sterics).
    """
    elec = np.linspace(1.0, 0.0, n_elec)
    vdw = np.linspace(1.0, 0.0, n_vdw)

    schedule: list[tuple[float, float]] = []
    for le in elec:
        schedule.append((float(le), 1.0))
    for ls in vdw[1:]:
        schedule.append((0.0, float(ls)))
    return schedule


# -----------------------------
# Reference relax (stability)
# -----------------------------
def relax_reference_system(
    system: mm.System,
    positions,
    platform_name: str = "CPU",
    temperature_k: float = 300.0,
    friction_per_ps: float = 10.0,
    step_fs: float = 0.5,
    minimize_max_iter: int = 5000,
    nsteps_warm: int = 2000,
) -> tuple:
    """Minimize + short warmup on the NON-alchemical reference system."""
    temperature = temperature_k * unit.kelvin
    integrator = mm.LangevinMiddleIntegrator(
        temperature,
        friction_per_ps / unit.picosecond,
        step_fs * unit.femtosecond,
    )
    integrator.setConstraintTolerance(1e-6)

    platform = mm.Platform.getPlatformByName(platform_name)
    context = mm.Context(system, integrator, platform)
    context.setPositions(positions)

    mm.LocalEnergyMinimizer.minimize(context, maxIterations=minimize_max_iter)

    st = context.getState(getEnergy=True, getPositions=True)
    E = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    pos_nm = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    print(
        f"DEBUG ref post-minimize: E(kJ/mol)={E:.3f}  pos_finite={np.isfinite(pos_nm).all()}  "
        f"max|pos|={np.max(np.abs(pos_nm)):.3f} nm"
    )
    if not np.isfinite(pos_nm).all():
        raise RuntimeError("Reference system has non-finite positions after minimization.")

    context.setVelocitiesToTemperature(temperature)
    integrator.step(nsteps_warm)

    st2 = context.getState(getPositions=True, getVelocities=True)
    return st2.getPositions(), st2.getVelocities()


# -------------------------------------------
# MBAR sampling (1 sampling ctx + 1 eval ctx)
# -------------------------------------------
def run_lambda_windows_mbar(
    alch_system: mm.System,
    positions,
    velocities,
    schedule: list[tuple[float, float]],
    platform_name: str = "CPU",
    temperature_k: float = 300.0,
    friction_per_ps: float = 20.0,
    step_fs: float = 0.25,
    nsteps_eq_each: int = 1500,
    nsteps_prod_each: int = 3000,
    sample_interval: int = 200,
    random_seed: int = 2025,
):
    """
    Stable MBAR collection:
      - sampling context integrates at state k
      - eval context re-evaluates energies at all states l for each sampled frame
    """
    temperature = temperature_k * unit.kelvin
    kT = (unit.MOLAR_GAS_CONSTANT_R * temperature).value_in_unit(unit.kilojoule_per_mole)
    beta = 1.0 / kT

    platform = mm.Platform.getPlatformByName(platform_name)
    lambda_names = _lambda_param_names(alch_system)
    print("DEBUG lambda params:", lambda_names)

    if "lambda_electrostatics" not in lambda_names or "lambda_sterics" not in lambda_names:
        raise RuntimeError(f"Expected lambda_electrostatics & lambda_sterics, got: {lambda_names}")

    integ_s = mm.LangevinMiddleIntegrator(
        temperature, friction_per_ps / unit.picosecond, step_fs * unit.femtosecond
    )
    integ_s.setConstraintTolerance(1e-6)
    integ_s.setRandomNumberSeed(int(random_seed))

    ctx_s = mm.Context(alch_system, integ_s, platform)
    ctx_s.setPositions(positions)
    if velocities is not None:
        ctx_s.setVelocities(velocities)
    else:
        ctx_s.setVelocitiesToTemperature(temperature)

    integ_e = mm.VerletIntegrator(1.0 * unit.femtosecond)
    ctx_e = mm.Context(alch_system, integ_e, platform)
    ctx_e.setPositions(positions)

    def set_state(ctx: mm.Context, le: float, ls: float):
        ctx.setParameter("lambda_electrostatics", float(le))
        ctx.setParameter("lambda_sterics", float(ls))

    K = len(schedule)
    n_samples = max(1, nsteps_prod_each // sample_interval)
    u_kln = np.zeros((K, K, n_samples), dtype=float)

    le0, ls0 = schedule[0]
    set_state(ctx_s, le0, ls0)
    mm.LocalEnergyMinimizer.minimize(ctx_s, maxIterations=2000)

    for k, (le_k, ls_k) in enumerate(schedule):
        print(f"DEBUG equil le={le_k:.3f} ls={ls_k:.3f}")
        set_state(ctx_s, le_k, ls_k)
        integ_s.step(nsteps_eq_each)

        print(f"DEBUG production le={le_k:.3f} ls={ls_k:.3f}  samples={n_samples}")
        for n in range(n_samples):
            integ_s.step(sample_interval)
            pos = ctx_s.getState(getPositions=True).getPositions()

            ctx_e.setPositions(pos)
            for l, (le_l, ls_l) in enumerate(schedule):
                set_state(ctx_e, le_l, ls_l)
                E = ctx_e.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                red_u = beta * E
                if not np.isfinite(red_u):
                    raise RuntimeError(
                        f"NaN reduced potential at k={k}, l={l}, "
                        f"(le_k,ls_k)=({le_k},{ls_k}), (le_l,ls_l)=({le_l},{ls_l})"
                    )
                u_kln[k, l, n] = red_u

    N_k = np.array([n_samples] * K, dtype=int)
    return u_kln, N_k, kT


def compute_deltaG_from_mbar(u_kln, N_k, kT_kJ_per_mol: float):
    """Compute ΔG between first and last state in schedule using MBAR."""
    mbar = MBAR(u_kln, N_k, verbose=False)
    res = mbar.compute_free_energy_differences()
    Delta_f = res["Delta_f"]
    dDelta_f = res["dDelta_f"]

    i = 0
    j = u_kln.shape[0] - 1
    dG = float(Delta_f[i, j] * kT_kJ_per_mol)
    dG_err = float(dDelta_f[i, j] * kT_kJ_per_mol)
    return dG, dG_err


def save_mbar_qc(u_kln, N_k, kT_kJ_per_mol: float, schedule, out_prefix: str):
    """
    Save minimal QC from MBAR:
      - Adjacent ΔG_i->i+1 and uncertainties
      - Cumulative ΔG from state0 to state i
    Outputs:
      {out_prefix}_mbar_qc.csv
      {out_prefix}_mbar_qc.png
      {out_prefix}_mbar_qc_cum.png
    """
    import matplotlib.pyplot as plt

    mbar = MBAR(u_kln, N_k, verbose=False)
    res = mbar.compute_free_energy_differences()
    Delta_f = res["Delta_f"]
    dDelta_f = res["dDelta_f"]

    K = u_kln.shape[0]

    dG_adj = np.array([Delta_f[i, i + 1] * kT_kJ_per_mol for i in range(K - 1)], dtype=float)
    dG_adj_err = np.array([dDelta_f[i, i + 1] * kT_kJ_per_mol for i in range(K - 1)], dtype=float)

    dG_cum = np.array([Delta_f[0, i] * kT_kJ_per_mol for i in range(K)], dtype=float)
    dG_cum_err = np.array([dDelta_f[0, i] * kT_kJ_per_mol for i in range(K)], dtype=float)

    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    csv_path = f"{out_prefix}_mbar_qc.csv"

    header = [
        "i",
        "le_i", "ls_i",
        "le_ip1", "ls_ip1",
        "dG_adj_kJmol", "dG_adj_err_kJmol",
        "dG_cum_to_i_kJmol", "dG_cum_err_to_i_kJmol",
    ]

    with open(csv_path, "w") as f:
        f.write(",".join(header) + "\n")
        for i in range(K):
            le_i, ls_i = schedule[i]
            if i < K - 1:
                le_j, ls_j = schedule[i + 1]
                row = [
                    i,
                    f"{le_i:.6f}", f"{ls_i:.6f}",
                    f"{le_j:.6f}", f"{ls_j:.6f}",
                    f"{dG_adj[i]:.6f}", f"{dG_adj_err[i]:.6f}",
                    f"{dG_cum[i]:.6f}", f"{dG_cum_err[i]:.6f}",
                ]
            else:
                # last state has no adjacent (i->i+1)
                row = [
                    i,
                    f"{le_i:.6f}", f"{ls_i:.6f}",
                    "", "",
                    "", "",
                    f"{dG_cum[i]:.6f}", f"{dG_cum_err[i]:.6f}",
                ]
            f.write(",".join(map(str, row)) + "\n")

    # plots
    png_path = f"{out_prefix}_mbar_qc.png"
    x_adj = np.arange(K - 1)
    x_cum = np.arange(K)

    plt.figure()
    plt.errorbar(x_adj, dG_adj, yerr=dG_adj_err, fmt="o")
    plt.xlabel("Adjacent window index i (i → i+1)")
    plt.ylabel("ΔG_adj (kJ/mol)")
    plt.title("MBAR adjacent window contributions")
    plt.tight_layout()
    plt.savefig(png_path, dpi=200)
    plt.close()

    cum_path = f"{out_prefix}_mbar_qc_cum.png"
    plt.figure()
    plt.errorbar(x_cum, dG_cum, yerr=dG_cum_err, fmt="o")
    plt.xlabel("Window index i (state 0 → i)")
    plt.ylabel("Cumulative ΔG (kJ/mol)")
    plt.title("MBAR cumulative free energy")
    plt.tight_layout()
    plt.savefig(cum_path, dpi=200)
    plt.close()

    print(f"Wrote QC CSV:   {csv_path}")
    print(f"Wrote QC plots: {png_path} and {cum_path}")


# -----------------------------
#        Run one ligand
# -----------------------------
def run_one_ligand(
    ligand_sdf: str,
    outdir: str,
    platform: str = "CPU",
    schedule: list[tuple[float, float]] | None = None,
    random_seed: int = 2025,
    do_qc: bool = True,
    # runtime knobs
    friction_per_ps: float = 20.0,
    step_fs: float = 0.25,
    nsteps_eq_each: int = 1500,
    nsteps_prod_each: int = 3000,
    sample_interval: int = 200,
):
    if schedule is None:
        schedule = make_decouple_schedule(n_elec=6, n_vdw=12)

    # make a tmp 3D sdf under outdir to avoid collisions
    os.makedirs(outdir, exist_ok=True)
    lig_base = os.path.splitext(os.path.basename(ligand_sdf))[0]
    tmp_sdf = os.path.join(outdir, f"{lig_base}.3d.sdf")
    ref_system, topology, positions = build_solvated_system(ligand_sdf, tmp_sdf=tmp_sdf)

    lig_atoms = get_ligand_atom_indices(topology)

    pos_relaxed, vel_relaxed = relax_reference_system(ref_system, positions, platform_name=platform)

    alch_sys = alchemical_system(ref_system, lig_atoms)

    u_kln, N_k, kT = run_lambda_windows_mbar(
        alch_sys,
        pos_relaxed,
        vel_relaxed,
        schedule=schedule,
        platform_name=platform,
        random_seed=random_seed,
        friction_per_ps=friction_per_ps,
        step_fs=step_fs,
        nsteps_eq_each=nsteps_eq_each,
        nsteps_prod_each=nsteps_prod_each,
        sample_interval=sample_interval,
    )

    dG_kJ, dG_err_kJ = compute_deltaG_from_mbar(u_kln, N_k, kT)

    outdir_seed = os.path.join(outdir, f"seed{random_seed}")
    os.makedirs(outdir_seed, exist_ok=True)

    result = {
        "ligand_sdf": ligand_sdf,
        "dG_hyd_kJmol": float(dG_kJ),
        "dG_sem_kJmol": float(dG_err_kJ),
        "schedule": schedule,
        "n_samples_per_state": int(N_k[0]),
        "random_seed": int(random_seed),
        "platform": platform,
        "temperature_K": 300.0,
        "friction_per_ps": float(friction_per_ps),
        "step_fs": float(step_fs),
        "nsteps_eq_each": int(nsteps_eq_each),
        "nsteps_prod_each": int(nsteps_prod_each),
        "sample_interval": int(sample_interval),
    }

    with open(os.path.join(outdir_seed, "hydration_result.json"), "w") as f:
        json.dump(result, f, indent=2)

    if do_qc:
        qc_dir = os.path.join("results", "qc")
        os.makedirs(qc_dir, exist_ok=True)
        lig_tag = os.path.basename(outdir.rstrip("/"))
        qc_prefix = os.path.join(qc_dir, f"{lig_tag}_seed{random_seed}")
        save_mbar_qc(u_kln, N_k, kT, schedule, qc_prefix)

    return dG_kJ, dG_err_kJ


def main():
    ligA = "data/ligands/BNZ.sdf"
    ligB = "data/ligands/MBN.sdf"

    outA = "results/hydration_abs/BNZ"
    outB = "results/hydration_abs/MBN"

    seeds = [2025, 2026, 2027]

    rows = []
    for seed in seeds:
        print(f"\n===== SEED {seed} =====")

        print("Running absolute hydration for ligA ...")
        dG_A, sem_A = run_one_ligand(ligA, outA, platform="CPU", random_seed=seed)

        print("Running absolute hydration for ligB ...")
        dG_B, sem_B = run_one_ligand(ligB, outB, platform="CPU", random_seed=seed)

        ddG = dG_B - dG_A
        ddG_sem = float(np.sqrt(sem_A**2 + sem_B**2))
        rows.append((seed, dG_A, sem_A, dG_B, sem_B, ddG, ddG_sem))

        print(f"ΔG_hyd(A) = {dG_A:.3f} ± {sem_A:.3f} kJ/mol")
        print(f"ΔG_hyd(B) = {dG_B:.3f} ± {sem_B:.3f} kJ/mol")
        print(f"ΔΔG_hyd(B-A) = {ddG:.3f} ± {ddG_sem:.3f} kJ/mol")

    ddGs = np.array([r[5] for r in rows], dtype=float)
    mean = float(ddGs.mean())
    std = float(ddGs.std(ddof=1)) if len(ddGs) > 1 else 0.0
    sem = float(std / np.sqrt(len(ddGs))) if len(ddGs) > 1 else 0.0

    os.makedirs("results", exist_ok=True)
    out_csv = "results/ddg_summary_repeats.csv"
    with open(out_csv, "w") as f:
        f.write("seed,dG_A_kJmol,sem_A_kJmol,dG_B_kJmol,sem_B_kJmol,ddG_BminusA_kJmol,ddG_sem_kJmol\n")
        for (seed, dG_A, sem_A, dG_B, sem_B, ddG, ddG_sem) in rows:
            f.write(f"{seed},{dG_A:.3f},{sem_A:.3f},{dG_B:.3f},{sem_B:.3f},{ddG:.3f},{ddG_sem:.3f}\n")
        f.write(f"MEAN,,,,,{mean:.3f},\n")
        f.write(f"STD,,,,,{std:.3f},\n")
        f.write(f"SEM,,,,,{sem:.3f},\n")

    print("\n===== REPEATS SUMMARY =====")
    print(f"ΔΔG mean = {mean:.3f} kJ/mol")
    print(f"ΔΔG std  = {std:.3f} kJ/mol")
    print(f"ΔΔG SEM  = {sem:.3f} kJ/mol")
    print(f"Wrote: {out_csv}")


if __name__ == "__main__":
    main()
