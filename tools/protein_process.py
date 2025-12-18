"""
Author: Yonglan Liu
Created: 2025-12
Project: OpenMM-based Free Energy Pipeline (Absolute Hydration)
"""
"""
This utility cleans and standardizes a protein–ligand PDB for downstream MD/FEP workflows:
* Fixes common PDB issues with PDBFixer (nonstandard residues, missing atoms, hydrogens).
* Keeps only the specified ligand by residue name (--keep_resname).
* Removes all other heterogens robustly (waters/ions optional).
* Optionally restricts ligand selection by chain and/or residue id to avoid keeping the wrong copy when multiple ligands exist.
"""
"""
Usage example:
python tools/fix_pdb_keep_ligand.py \
        -i data/pdb/[PDBid].pdb \
        -o data/protein/[PDBid]_fixed_keep[ligand Resname].pdb \
        --keep_resname [ligand Resname] \
  
Choose the correct ligand when multiple copies exist:
python tools/fix_pdb_keep_ligand.py \
        -i data/pdb/[PDBid].pdb \
        -o data/protein/[PDBid]_fixed_keep[ligand Resname].pdb \
        --keep_resname [ligand Resname] \
        --keep_chain [chain id] \
        --keep_resid [residue id]

Keep waters and/or common ions:
By default, crystallographic waters and common ions are removed. Use flags to keep them:
# Keep waters:
python tools/fix_pdb_keep_ligand.py \
      -i data/pdb/[PDBid].pdb \
      -o data/protein/[PDBid]_fixed_keep[ligand Resname].pdb \
      --keep_resname [ligand Resname] \
      --keep_water

# Keep ions:
python tools/fix_pdb_keep_ligand.py \
      -i data/pdb/[PDBid].pdb \
      -o data/protein/[PDBid]_fixed_keep[ligand Resname].pdb \
      --keep_resname [ligand Resname] \
      --keep_ions

# Keep both waters and ions:
python tools/fix_pdb_keep_ligand.py \
      -i data/pdb/[PDBid].pdb \
      -o data/protein/[PDBid]_fixed_keep[ligand Resname].pdb \
        --keep_resname [ligand Resname] \
        --keep_water \
        --keep_ions
"""
#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile, Modeller

WATER_NAMES = {"HOH", "WAT", "SOL", "TIP3"}
COMMON_IONS = {"NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "CU", "CO", "NI"}

# 标准蛋白残基 + 常见质子化/端基/修饰名（MVP 足够稳）
STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET",
    "PHE","PRO","SER","THR","TRP","TYR","VAL",
    # Histidine tautomers / charged
    "HID","HIE","HIP",
    # Asp/Glu protonated
    "ASH","GLH",
    # Cys variants in some files
    "CYX","CYM",
    # Lys variants
    "LYN",
    # N/C terminal caps sometimes appear
    "ACE","NME"
}

# if there are DNA/RNA, you can add them here (default is not to keep, to avoid accidental inclusion of unwanted components)
STANDARD_NA = {
    "DA","DT","DG","DC","A","U","G","C","I"
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Fix protein PDB and keep only selected ligand by resname; remove all other heterogens robustly."
    )
    p.add_argument("-i", "--in_pdb", required=True, help="Input PDB path")
    p.add_argument("-o", "--out_pdb", required=True, help="Output PDB path")

    p.add_argument("--keep_resname", required=True, help="Ligand residue name to KEEP (e.g., BNZ, MBN, LIG)")
    p.add_argument("--keep_chain", default=None, help="Optional: keep ligand only on this chain (e.g., A)")
    p.add_argument("--keep_resid", default=None, help="Optional: keep ligand only with this residue id (e.g., 401)")

    p.add_argument("--ph", type=float, default=7.0, help="pH for adding hydrogens (default: 7.0)")

    p.add_argument("--keep_water", action="store_true", help="Keep crystallographic waters (default: remove waters)")
    p.add_argument("--keep_ions", action="store_true", help="Keep common ions (default: remove ions)")

    p.add_argument("--keep_nucleic_acids", action="store_true",
                   help="Keep nucleic acid residues (DA/DT/...); default: delete them as heterogens.")
    p.add_argument("--keep_all_matching_ligands", action="store_true",
                   help="Keep ALL ligands with keep_resname (default: keep only one unless chain/resid specified).")

    return p.parse_args()


def _rid(res) -> str:
    return str(res.id)


def main():
    args = parse_args()
    keep_resname = args.keep_resname.upper()

    # -----------------------
    # 1) PDBFixer: 修蛋白 + 加氢
    # -----------------------
    fixer = PDBFixer(filename=args.in_pdb)

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    fixer.findMissingResidues()
    chain_min = {}
    chain_max = {}
    for chain in fixer.topology.chains():
        ids = []
        for res in chain.residues():
            try:
                ids.append(int(res.id))  # PDB residue number
            except Exception:
                pass
        if ids:
            chain_min[chain.index] = min(ids)
            chain_max[chain.index] = max(ids)
    # PDBFixer 的 missingResidues: key=(chainIndex, resId) value=[(resName, resId), ...] 之类
    filtered = {}
    for key, missing_list in fixer.missingResidues.items():
        chain_idx, anchor_resid = key
        # if the missing segment is near terminal, consider it terminal missing and drop it
        if chain_idx in chain_min and chain_idx in chain_max:
            if anchor_resid <= chain_min[chain_idx] or anchor_resid >= chain_max[chain_idx]:
                continue  # drop terminal missing residues
        filtered[key] = missing_list
    fixer.missingResidues = filtered
       
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(pH=args.ph)

    modeller = Modeller(fixer.topology, fixer.positions)

    # -------------------------------------------------------------------
    # 2) search ligand (there might be several ligands with same resname)
    # -------------------------------------------------------------------
    lig_matches = [res for res in modeller.topology.residues() if res.name.upper() == keep_resname]
    if len(lig_matches) == 0:
        raise RuntimeError(
            f"No ligand residue with resname '{keep_resname}' found. "
            f"Please check HETATM lines in the input PDB."
        )

    # Only keep only one ligand to avoid FEP alchemizing two by mistake
    def match_filter(res):
        if args.keep_chain is not None and res.chain.id != args.keep_chain:
            return False
        if args.keep_resid is not None and _rid(res) != str(args.keep_resid):
            return False
        return True

    if args.keep_chain is not None or args.keep_resid is not None:
        selected = [r for r in lig_matches if match_filter(r)]
        if len(selected) == 0:
            raise RuntimeError(
                f"Found ligands {[(r.name, r.chain.id, _rid(r)) for r in lig_matches]} "
                f"but none match keep_chain={args.keep_chain} keep_resid={args.keep_resid}."
            )
        keep_ligs = selected if args.keep_all_matching_ligands else [selected[0]]
    else:
        keep_ligs = lig_matches if args.keep_all_matching_ligands else [lig_matches[0]]

    keep_lig_set = set(keep_ligs)

    # -----------------------
    # 3) construct delete list: delete water, ions, and all "non-standard protein/nucleic acid/kept ligand"
    # -----------------------
    delete_res = []
    for res in modeller.topology.residues():
        name = res.name.upper()

        # keep the target ligand (specified one or more)
        if name == keep_resname and res in keep_lig_set:
            continue

        # delete other ligands with same resname
        if name == keep_resname and res not in keep_lig_set:
            delete_res.append(res)
            continue

        # water
        if name in WATER_NAMES:
            if not args.keep_water:
                delete_res.append(res)
            continue

        # ions
        if name in COMMON_IONS:
            if not args.keep_ions:
                delete_res.append(res)
            continue

        # standard protein residues: keep
        if name in STANDARD_AA:
            continue

        # nucleic acids: keep if requested
        if args.keep_nucleic_acids and name in STANDARD_NA:
            continue

        # delete other heterogens (this is the key step for "clean" PDB)
        delete_res.append(res)

    if delete_res:
        modeller.delete(delete_res)

    # -----------------------
    # 4) write output PDB
    # -----------------------
    with open(args.out_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)

    # --------------------------------------------------------------
    # 5) print summary (you can use it to confirm deletion is clean)
    # --------------------------------------------------------------
    kept_after = [(r.name, r.chain.id, _rid(r)) for r in modeller.topology.residues() if r.name.upper() == keep_resname]
    other_het_after = []
    for r in modeller.topology.residues():
        n = r.name.upper()
        if n in STANDARD_AA or n in WATER_NAMES or n in COMMON_IONS or n == keep_resname:
            continue
        if args.keep_nucleic_acids and n in STANDARD_NA:
            continue
        other_het_after.append((r.name, r.chain.id, _rid(r)))

    print("Input :", args.in_pdb)
    print("Output:", args.out_pdb)
    print("Keep ligand resname:", keep_resname)
    print("Ligands found BEFORE:", [(r.name, r.chain.id, _rid(r)) for r in lig_matches])
    print("Ligands kept AFTER :", kept_after)
    print("Other heterogens AFTER (should be empty):", other_het_after)


if __name__ == "__main__":
    main()
