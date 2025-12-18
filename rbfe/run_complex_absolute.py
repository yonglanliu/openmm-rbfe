"""
run_complex_absolute.py
======================

Absolute complex-leg decoupling (protein-ligand complex) using OpenMM + openmmtools alchemy + MBAR.

Key stability features (important for real PDBs):
- Pre-minimize coordinate sanity check (finite positions)
- Damped dynamics before minimization (strong friction + small dt)
- Automatic fallback relax: if reference minimize hits NaN, rebuild a "looser"
  system (constraints=None) to relieve clashes, then transfer relaxed positions back
  to the normal system.

Outputs:
- results/complex_abs/<lig_tag>/seed<seed>/complex_result.json
- results/qc/<lig_tag>_complex_seed<seed>_mbar_qc.csv
- results/qc/<lig_tag>_complex_seed<seed>_mbar_qc.png
- results/qc/<lig_tag>_complex_seed<seed>_mbar_qc_cum.png
"""

import os
import json
import numpy as np

from openmm import unit
import openmm as mm
import openmm.app as app

from openff.toolkit import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator

from openmmtools.alchemy import AlchemicalRegion
from pymbar import MBAR


# -----------------------------
# Utilities: lambda parameters
# -----------------------------
def _lambda_param_names(system: mm.System) -> list[str]:
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
    Two-stage schedule:
      Stage 1: electrostatics 1 -> 0 (sterics = 1)
      Stage 2: sterics       1 -> 0 (electrostatics = 0)
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
#    Complex system build
# -----------------------------
def load_complex_from_pdb(complex_pdb: str) -> tuple[app.Topology, unit.Quantity]:
    pdb = app.PDBFile(complex_pdb)
    return pdb.topology, pdb.positions


def identify_ligand_atoms_by_resname(topology: app.Topology, ligand_resname: str) -> list[int]:
    target = ligand_resname.upper().strip()
    lig_atoms = []
    for atom in topology.atoms():
        if atom.residue.name.upper().strip() == target:
            lig_atoms.append(atom.index)
    if len(lig_atoms) == 0:
        raise RuntimeError(f"Could not find ligand atoms for resname='{ligand_resname}'.")
    return lig_atoms


def build_solvated_complex_system(
    complex_pdb: str,
    ligand_sdf: str,
    ligand_resname: str,
    padding_nm: float = 1.2,
    constraints_mode="HBonds",  # "HBonds" (default), "AllBonds", or None
):
    """
    Build solvated protein-ligand complex using AMBER protein FF + TIP3P + GAFF ligand.
    Returns: (system, topology, positions, lig_atom_indices)
    """
    topology, positions = load_complex_from_pdb(complex_pdb)

    offmol = Molecule.from_file(ligand_sdf)

    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    gaff = GAFFTemplateGenerator(molecules=[offmol], forcefield="gaff-2.11")
    forcefield.registerTemplateGenerator(gaff.generator)

    modeller = app.Modeller(topology, positions)
    modeller.addSolvent(forcefield, model="tip3p", padding=padding_nm * unit.nanometer)

    if constraints_mode == "HBonds":
        constraints = app.HBonds
    elif constraints_mode == "AllBonds":
        constraints = app.AllBonds
    elif constraints_mode is None or constraints_mode == "None":
        constraints = None
    else:
        raise ValueError(f"Unknown constraints_mode={constraints_mode}")

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=constraints,
        rigidWater=True,
    )

    lig_atoms = identify_ligand_atoms_by_resname(modeller.topology, ligand_resname)
    return system, modeller.topology, modeller.positions, lig_atoms


# -----------------------------
#     Restraints (distance)
# -----------------------------
def add_distance_restraint(
    system: mm.System,
    topology: app.Topology,
    positions,
    ligand_resname: str,
    k_kJ_per_mol_nm2: float = 2000.0,
    r0_nm: float | None = None,
) -> dict:
    """
    Harmonic distance restraint between nearest protein heavy atom and first ligand heavy atom.
    """
    lig_res = ligand_resname.upper().strip()
    water_names = {"HOH", "WAT", "SOL", "TIP3"}

    # Convert positions to numpy (N,3) in nm
    pos_nm = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)

    # pick ligand anchor: first heavy atom in ligand residue
    lig_heavy = []
    for a in topology.atoms():
        if a.residue.name.upper().strip() == lig_res and a.element is not None and a.element.atomic_number > 1:
            lig_heavy.append(a.index)
    if not lig_heavy:
        raise RuntimeError("No ligand heavy atoms found (anchor selection failed).")
    lig_anchor = int(lig_heavy[0])
    lig_xyz = pos_nm[lig_anchor]

    # protein heavy atoms (exclude water + ligand)
    prot_heavy = []
    for a in topology.atoms():
        rn = a.residue.name.upper().strip()
        if rn in water_names:
            continue
        if rn == lig_res:
            continue
        if a.element is None or a.element.atomic_number <= 1:
            continue
        prot_heavy.append(int(a.index))
    if not prot_heavy:
        raise RuntimeError("No protein heavy atoms found (anchor selection failed).")

    prot_xyz = pos_nm[prot_heavy]
    d2 = np.sum((prot_xyz - lig_xyz[None, :]) ** 2, axis=1)
    prot_anchor = prot_heavy[int(np.argmin(d2))]

    if r0_nm is None:
        r0_nm = float(np.sqrt(np.min(d2)))

    force = mm.CustomBondForce("0.5*k*(r-r0)^2")
    force.addGlobalParameter("k", float(k_kJ_per_mol_nm2))
    force.addGlobalParameter("r0", float(r0_nm))
    force.addBond(int(prot_anchor), int(lig_anchor), [])
    system.addForce(force)

    return {
        "type": "harmonic_distance",
        "k_kJ_per_mol_nm2": float(k_kJ_per_mol_nm2),
        "r0_nm": float(r0_nm),
        "protein_anchor_atom_index": int(prot_anchor),
        "ligand_anchor_atom_index": int(lig_anchor),
    }


# -----------------------------
#      Alchemical system
# -----------------------------
def alchemical_system(reference_system: mm.System, alchemical_atoms: list[int]) -> mm.System:
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
# Reference relax (robust)
# -----------------------------
def _positions_to_numpy_nm(positions) -> np.ndarray:
    return np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)

# Get potential energy from context in kJ/mol
def _context_energy_kjmol(context: mm.Context) -> float:
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

# Robust reference system relax
def relax_reference_system_robust(
    system: mm.System,
    positions,
    platform_name: str = "CPU",
    temperature_k: float = 300.0,
    random_seed: int = 2025,
    minimize_max_iter: int = 5000,
    damp_steps: int = 100,
    damp_friction_per_ps: float = 100.0,
    damp_step_fs: float = 0.25,
    warm_steps: int = 500,
) -> tuple:
    """
    Robust relax:
      1) check finite positions
      2) damped dynamics to relieve clashes
      3) minimize
      4) short warmup
    Returns (positions, velocities)
    """
    # pre-check：if positions are finite orginally
    pos_nm = _positions_to_numpy_nm(positions)
    print("DEBUG pre-min pos finite:", np.isfinite(pos_nm).all(), "max|pos|:", float(np.max(np.abs(pos_nm))))
    if not np.isfinite(pos_nm).all():
        bad = np.where(~np.isfinite(pos_nm))
        raise RuntimeError(f"Non-finite positions before minimize! bad indices: {bad}")

    temperature = temperature_k * unit.kelvin
    platform = mm.Platform.getPlatformByName(platform_name)

    # use LangevinMiddle (stable)
    integrator = mm.LangevinMiddleIntegrator(
        temperature,
        damp_friction_per_ps / unit.picosecond,
        damp_step_fs * unit.femtosecond,
    )
    integrator.setConstraintTolerance(1e-6)
    integrator.setRandomNumberSeed(int(random_seed))

    context = mm.Context(system, integrator, platform)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature)

    # quick initial energy (may be huge, but should be finite)
    try:
        E0 = _context_energy_kjmol(context)
        print(f"DEBUG pre-min energy (kJ/mol): {E0:.3e}")
    except Exception as e:
        print("DEBUG pre-min energy eval failed:", repr(e))

    # damped dynamics before minimize
    if damp_steps > 0:
        integrator.step(int(damp_steps))

    # minimize
    mm.LocalEnergyMinimizer.minimize(context, maxIterations=int(minimize_max_iter))

    st = context.getState(getEnergy=True, getPositions=True)
    E = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    pos_nm2 = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    print(
        f"DEBUG post-min: E(kJ/mol)={E:.3f}  pos_finite={np.isfinite(pos_nm2).all()}  max|pos|={np.max(np.abs(pos_nm2)):.3f} nm"
    )
    if not np.isfinite(pos_nm2).all():
        raise RuntimeError("Non-finite positions after minimization.")

    # warmup
    integrator.step(int(warm_steps))
    st2 = context.getState(getPositions=True, getVelocities=True)
    return st2.getPositions(), st2.getVelocities()


def relax_reference_with_fallback(
    build_fn,
    build_kwargs: dict,
    platform_name: str,
    temperature_k: float,
    random_seed: int,
) -> tuple:
    """
    Try relax on normal system; if NaN during minimize, rebuild a looser system (constraints=None),
    relax there, then return relaxed positions/velocities from the looser relax. The caller can then
    seed the normal system contexts with these relaxed positions.
    """
    # 1) normal relax (HBonds is default for build in this script)
    system, topology, positions, lig_atoms, restraint_info = build_fn(**build_kwargs)
    try:
        pos_relaxed, vel_relaxed = relax_reference_system_robust(
            system=system,
            positions=positions,
            platform_name=platform_name,
            temperature_k=temperature_k,
            random_seed=random_seed,
        )
        return system, topology, positions, lig_atoms, restraint_info, pos_relaxed, vel_relaxed
    except Exception as e:
        msg = str(e)
        if "Particle coordinate is NaN" not in msg and "Non-finite" not in msg:
            raise  # not the NaN class we are handling
        print("\nWARNING: reference relax failed with NaN-like error.")
        print("         Trying fallback relax with constraints=None (looser system) ...\n")

    # 2) fallback: constraints=None
    build_kwargs2 = dict(build_kwargs)
    build_kwargs2["constraints_mode"] = None
    system2, topology2, positions2, lig_atoms2, restraint_info2 = build_fn(**build_kwargs2)

    pos_relaxed2, vel_relaxed2 = relax_reference_system_robust(
        system=system2,
        positions=positions2,
        platform_name=platform_name,
        temperature_k=temperature_k,
        random_seed=random_seed,
        damp_steps=200,
        damp_friction_per_ps=200.0,
        damp_step_fs=0.2,
        minimize_max_iter=8000,
        warm_steps=1000,
    )

    print("INFO: fallback relax succeeded (constraints=None).")
    # return the *normal* system for production but the fallback relaxed positions/vels
    # Caller will build normal system again and seed with relaxed coords.
    return None, None, None, None, None, pos_relaxed2, vel_relaxed2


# -------------------------------------------
# MBAR sampling (sampling ctx + eval ctx)
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

    # minimize once at initial state
    le0, ls0 = schedule[0]
    set_state(ctx_s, le0, ls0)
    mm.LocalEnergyMinimizer.minimize(ctx_s, maxIterations=2000)

    for k, (le_k, ls_k) in enumerate(schedule):
        print(f"DEBUG equil le={le_k:.3f} ls={ls_k:.3f}")
        set_state(ctx_s, le_k, ls_k)
        ctx_s.setVelocitiesToTemperature(temperature)  # re-thermalize per state helps stability
        integ_s.step(int(nsteps_eq_each))

        print(f"DEBUG production le={le_k:.3f} ls={ls_k:.3f}  samples={n_samples}")
        for n in range(n_samples):
            integ_s.step(int(sample_interval))
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

    with open(csv_path, "w") as f:
        f.write("i,le_i,ls_i,le_ip1,ls_ip1,dG_adj_kJmol,dG_adj_err_kJmol,dG_cum_to_i_kJmol,dG_cum_err_to_i_kJmol\n")
        for i in range(K):
            le_i, ls_i = schedule[i]
            if i < K - 1:
                le_j, ls_j = schedule[i + 1]
                f.write(
                    f"{i},{le_i:.6f},{ls_i:.6f},{le_j:.6f},{ls_j:.6f},"
                    f"{dG_adj[i]:.6f},{dG_adj_err[i]:.6f},{dG_cum[i]:.6f},{dG_cum_err[i]:.6f}\n"
                )
            else:
                f.write(f"{i},{le_i:.6f},{ls_i:.6f},,,,{dG_cum[i]:.6f},{dG_cum_err[i]:.6f}\n")

    png_path = f"{out_prefix}_mbar_qc.png"
    x_adj = np.arange(K - 1)

    plt.figure()
    plt.errorbar(x_adj, dG_adj, yerr=dG_adj_err, fmt="o")
    plt.xlabel("Adjacent window index i (i → i+1)")
    plt.ylabel("ΔG_adj (kJ/mol)")
    plt.title("MBAR adjacent window contributions (complex leg)")
    plt.tight_layout()
    plt.savefig(png_path, dpi=200)
    plt.close()

    cum_path = f"{out_prefix}_mbar_qc_cum.png"
    x_cum = np.arange(K)

    plt.figure()
    plt.errorbar(x_cum, dG_cum, yerr=dG_cum_err, fmt="o")
    plt.xlabel("Window index i (state 0 → i)")
    plt.ylabel("Cumulative ΔG (kJ/mol)")
    plt.title("MBAR cumulative free energy (complex leg)")
    plt.tight_layout()
    plt.savefig(cum_path, dpi=200)
    plt.close()

    print(f"Wrote QC CSV:   {csv_path}")
    print(f"Wrote QC plots: {png_path} and {cum_path}")


# -----------------------------
# Build wrapper (so fallback can rebuild)
# -----------------------------
def _build_reference_complex(
    complex_pdb: str,
    ligand_sdf: str,
    ligand_resname: str,
    padding_nm: float,
    constraints_mode,
    restraint_k_kJ_per_mol_nm2: float,
    restraint_r0_nm: float | None,
):
    system, topology, positions, lig_atoms = build_solvated_complex_system(
        complex_pdb=complex_pdb,
        ligand_sdf=ligand_sdf,
        ligand_resname=ligand_resname,
        padding_nm=padding_nm,
        constraints_mode=constraints_mode,
    )
    restraint_info = add_distance_restraint(
        system=system,
        topology=topology,
        positions=positions,
        ligand_resname=ligand_resname,
        k_kJ_per_mol_nm2=restraint_k_kJ_per_mol_nm2,
        r0_nm=restraint_r0_nm,
    )
    return system, topology, positions, lig_atoms, restraint_info


# -----------------------------
# Run one complex
# -----------------------------
def run_one_complex(
    complex_pdb: str,
    ligand_sdf: str,
    ligand_resname: str,
    outdir: str,
    platform: str = "CPU",
    schedule: list[tuple[float, float]] | None = None,
    random_seed: int = 2025,
    do_qc: bool = True,
    # runtime knobs
    temperature_k: float = 300.0,
    friction_per_ps: float = 20.0,
    step_fs: float = 0.25,
    nsteps_eq_each: int = 1500,
    nsteps_prod_each: int = 3000,
    sample_interval: int = 200,
    # build knobs
    padding_nm: float = 1.2,
    constraints_mode="HBonds",
    # restraint knobs
    restraint_k_kJ_per_mol_nm2: float = 2000.0,
    restraint_r0_nm: float | None = None,
):
    if schedule is None:
        schedule = make_decouple_schedule(n_elec=6, n_vdw=12)

    os.makedirs(outdir, exist_ok=True)

    # 1) robust reference relax with fallback
    build_kwargs = dict(
        complex_pdb=complex_pdb,
        ligand_sdf=ligand_sdf,
        ligand_resname=ligand_resname,
        padding_nm=padding_nm,
        constraints_mode=constraints_mode,
        restraint_k_kJ_per_mol_nm2=restraint_k_kJ_per_mol_nm2,
        restraint_r0_nm=restraint_r0_nm,
    )

    sys1, top1, pos1, lig_atoms1, restr1, pos_relaxed, vel_relaxed = relax_reference_with_fallback(
        build_fn=_build_reference_complex,
        build_kwargs=build_kwargs,
        platform_name=platform,
        temperature_k=temperature_k,
        random_seed=random_seed,
    )

    # 2) build the "normal" reference system fresh (so we have the correct System object)
    ref_system, topology, positions, lig_atoms, restraint_info = _build_reference_complex(**build_kwargs)

    # 3) seed with relaxed positions from robust relax
    # (positions from fallback may not be bitwise-identical topology, but atom ordering is the same
    #  because it's from the same modeller build path; if not, we'll find out quickly in minimize below)
    positions_seed = pos_relaxed
    velocities_seed = vel_relaxed

    # 4) build alchemical system
    alch_sys = alchemical_system(ref_system, lig_atoms)

    # 5) run MBAR sampling
    u_kln, N_k, kT = run_lambda_windows_mbar(
        alch_system=alch_sys,
        positions=positions_seed,
        velocities=velocities_seed,
        schedule=schedule,
        platform_name=platform,
        temperature_k=temperature_k,
        friction_per_ps=friction_per_ps,
        step_fs=step_fs,
        nsteps_eq_each=nsteps_eq_each,
        nsteps_prod_each=nsteps_prod_each,
        sample_interval=sample_interval,
        random_seed=random_seed,
    )

    dG_kJ, dG_err_kJ = compute_deltaG_from_mbar(u_kln, N_k, kT)

    outdir_seed = os.path.join(outdir, f"seed{random_seed}")
    os.makedirs(outdir_seed, exist_ok=True)

    result = {
        "complex_pdb": complex_pdb,
        "ligand_sdf": ligand_sdf,
        "ligand_resname": ligand_resname,
        "dG_complex_kJmol": float(dG_kJ),
        "dG_sem_kJmol": float(dG_err_kJ),
        "schedule": schedule,
        "n_samples_per_state": int(N_k[0]),
        "random_seed": int(random_seed),
        "platform": platform,
        "temperature_K": float(temperature_k),
        "friction_per_ps": float(friction_per_ps),
        "step_fs": float(step_fs),
        "nsteps_eq_each": int(nsteps_eq_each),
        "nsteps_prod_each": int(nsteps_prod_each),
        "sample_interval": int(sample_interval),
        "constraints_mode": str(constraints_mode),
        "padding_nm": float(padding_nm),
        "restraint": restraint_info,
        "alchemical_atoms_count": int(len(lig_atoms)),
    }

    with open(os.path.join(outdir_seed, "complex_result.json"), "w") as f:
        json.dump(result, f, indent=2)

    if do_qc:
        qc_dir = os.path.join("results", "qc")
        os.makedirs(qc_dir, exist_ok=True)

        lig_tag = os.path.basename(outdir.rstrip("/"))
        qc_prefix = os.path.join(qc_dir, f"{lig_tag}_complex_seed{random_seed}")
        save_mbar_qc(u_kln, N_k, kT, schedule, qc_prefix)

    return dG_kJ, dG_err_kJ


def main():
    # ---- EDIT THESE 3 ----
    complex_pdb = "data/protein/4w53_fixed_keepMBN.pdb"
    ligand_sdf = "data/ligands/MBN.sdf"
    ligand_resname = "MBN"

    outdir = "results/complex_abs/MBN"
    seeds = [2025, 2026, 2027]
    schedule = make_decouple_schedule(n_elec=6, n_vdw=12)

    rows = []
    for seed in seeds:
        print(f"\n===== SEED {seed} =====")
        dG, sem = run_one_complex(
            complex_pdb=complex_pdb,
            ligand_sdf=ligand_sdf,
            ligand_resname=ligand_resname,
            outdir=outdir,
            platform="CPU",
            schedule=schedule,
            random_seed=seed,
            do_qc=True,
            # runtime knobs (stable defaults)
            temperature_k=300.0,
            friction_per_ps=20.0,
            step_fs=0.25,
            nsteps_eq_each=1500,
            nsteps_prod_each=3000,
            sample_interval=200,
            # build knobs
            padding_nm=1.2,
            constraints_mode="HBonds",
            # restraint knobs
            restraint_k_kJ_per_mol_nm2=2000.0,
            restraint_r0_nm=None,
        )
        print(f"ΔG_complex = {dG:.3f} ± {sem:.3f} kJ/mol")
        rows.append((seed, dG, sem))

    dGs = np.array([r[1] for r in rows], dtype=float)
    mean = float(dGs.mean())
    std = float(dGs.std(ddof=1)) if len(dGs) > 1 else 0.0
    sem = float(std / np.sqrt(len(dGs))) if len(dGs) > 1 else 0.0

    os.makedirs("results", exist_ok=True)
    out_csv = f"results/complex_summary_repeats.csv"
    with open(out_csv, "w") as f:
        f.write("seed,dG_complex_kJmol,sem_kJmol\n")
        for seed, dG, s in rows:
            f.write(f"{seed},{dG:.6f},{s:.6f}\n")
        f.write(f"MEAN,{mean:.6f},\n")
        f.write(f"STD,{std:.6f},\n")
        f.write(f"SEM,{sem:.6f},\n")

    print("\n===== REPEATS SUMMARY =====")
    print(f"ΔG mean = {mean:.3f} kJ/mol")
    print(f"ΔG std  = {std:.3f} kJ/mol")
    print(f"ΔG SEM  = {sem:.3f} kJ/mol")
    print(f"Wrote: {out_csv}")


if __name__ == "__main__":
    main()
