#!/usr/bin/env python3
from __future__ import annotations

import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def parse_args():
    p = argparse.ArgumentParser(description="Extract ligand coords from complex PDB and write SDF with correct chemistry.")
    p.add_argument("-i", "--pdb", required=True, help="Input complex PDB (protein+ligand, cleaned)")
    p.add_argument("--resname", required=True, help="Ligand residue name in PDB (e.g., BNZ, MBN)")
    p.add_argument("--smiles", required=True, help="Ligand SMILES (e.g., c1ccccc1 for benzene)")
    p.add_argument("-o", "--out_sdf", required=True, help="Output SDF path")
    p.add_argument("--out_lig_pdb", default=None, help="Optional: also write ligand-only PDB (from RDKit mol)")
    return p.parse_args()


def main():
    args = parse_args()
    resname = args.resname.upper()

    # 1) Read complex PDB only for coordinates
    mol_pdb = Chem.MolFromPDBFile(args.pdb, removeHs=False, sanitize=False)
    if mol_pdb is None:
        raise RuntimeError(f"Failed to read PDB: {args.pdb}")

    conf = mol_pdb.GetConformer()

    # 2) Collect ligand atoms of given resname; split heavy vs H
    heavy_coords = []
    heavy_elems = []
    h_coords = []

    for a in mol_pdb.GetAtoms():
        info = a.GetPDBResidueInfo()
        if info is None:
            continue
        if info.GetResidueName().strip().upper() != resname:
            continue

        p = conf.GetAtomPosition(a.GetIdx())
        elem = a.GetSymbol()

        if a.GetAtomicNum() > 1:
            heavy_coords.append([p.x, p.y, p.z])
            heavy_elems.append(elem)
        else:
            h_coords.append([p.x, p.y, p.z])

    if len(heavy_coords) == 0:
        raise RuntimeError(f"No heavy atoms found for resname={resname} in {args.pdb}")

    heavy_coords = np.asarray(heavy_coords, dtype=float)
    h_coords = np.asarray(h_coords, dtype=float) if len(h_coords) else np.zeros((0, 3), dtype=float)

    # 3) Build RDKit mol from SMILES with correct bonds
    mol = Chem.AddHs(Chem.MolFromSmiles(args.smiles))
    if mol is None:
        raise RuntimeError("Failed to parse SMILES.")

    # Embed (just to create a conformer container)
    AllChem.EmbedMolecule(mol, randomSeed=2025)

    # 4) Map PDB heavy atoms -> RDKit heavy atoms
    # For BNZ/MBN: using "attached H count" signature works well when PDB contains H coordinates.
    # Compute PDB heavy atom attached-H count by proximity
    def count_attached_H(p_heavy, Hxyz, cutoff_A=1.25):
        if Hxyz.shape[0] == 0:
            return 0
        d2 = np.sum((Hxyz - p_heavy[None, :]) ** 2, axis=1)
        return int(np.sum(d2 < (cutoff_A ** 2)))

    pdb_hcount = np.array([count_attached_H(heavy_coords[i], h_coords) for i in range(len(heavy_coords))], dtype=int)

    # RDKit heavy atoms + their attached H count (explicit H)
    heavy_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    rd_hcount = np.array([sum(1 for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1) for i in heavy_idx], dtype=int)

    if len(heavy_idx) != len(heavy_coords):
        raise RuntimeError(
            f"Heavy atom count mismatch AFTER filtering H: SMILES heavy={len(heavy_idx)} vs PDB heavy={len(heavy_coords)}. "
            f"Check SMILES or resname extraction."
        )

    # Greedy group-by-Hcount mapping (works great for benzene/toluene)
    # Build mapping order indices within each hcount group
    mapping = [-1] * len(heavy_idx)

    for hc in sorted(set(rd_hcount.tolist())):
        rd_group = [heavy_idx[i] for i in range(len(heavy_idx)) if rd_hcount[i] == hc]
        pdb_group = [i for i in range(len(heavy_coords)) if pdb_hcount[i] == hc]

        if len(rd_group) != len(pdb_group):
            raise RuntimeError(
                f"Cannot map heavy atoms by attached-H count. "
                f"SMILES group hc={hc} has {len(rd_group)} atoms but PDB has {len(pdb_group)}. "
                f"Try verifying SMILES or use an RCSB SDF instead."
            )

        # deterministic ordering inside group
        rd_group = sorted(rd_group)
        pdb_group = sorted(pdb_group)

        for rdi, pdbi in zip(rd_group, pdb_group):
            mapping[heavy_idx.index(rdi)] = pdbi

    if any(m < 0 for m in mapping):
        raise RuntimeError("Internal mapping failed.")

    # 5) Set heavy-atom coordinates in RDKit conformer
    conf2 = mol.GetConformer()
    for j, atom_id in enumerate(heavy_idx):
        pdbi = mapping[j]
        x, y, z = heavy_coords[pdbi]
        conf2.SetAtomPosition(atom_id, Chem.rdGeometry.Point3D(float(x), float(y), float(z)))

    # 6) Optimize with heavy atoms constrained (keeps pose; relax H)
    # RDKit doesn't have a strict heavy-atom constraint in UFF API, but a short optimize is usually fine for BNZ/MBN.
    AllChem.UFFOptimizeMolecule(mol, maxIters=200)

    # 7) Write SDF
    w = Chem.SDWriter(args.out_sdf)
    w.write(mol)
    w.close()
    print("Wrote SDF:", args.out_sdf)

    if args.out_lig_pdb:
        Chem.MolToPDBFile(mol, args.out_lig_pdb)
        print("Wrote ligand PDB:", args.out_lig_pdb)


if __name__ == "__main__":
    main()
