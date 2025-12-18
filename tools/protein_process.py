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

# 若你遇到 DNA/RNA，可把这类加进去（默认先不保留，以免误保留杂质）
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
        # 只要缺失段发生在链的两端附近，就认为是 terminal missing，丢掉
        # （anchor_resid 在已存在范围之外/紧贴边界，通常就是 terminal gap）
        if chain_idx in chain_min and chain_idx in chain_max:
            if anchor_resid <= chain_min[chain_idx] or anchor_resid >= chain_max[chain_idx]:
                continue  # drop terminal missing residues
        filtered[key] = missing_list
    fixer.missingResidues = filtered
       
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(pH=args.ph)

    modeller = Modeller(fixer.topology, fixer.positions)

    # -----------------------
    # 2) 找 ligand（同 resname 可能有多个）
    # -----------------------
    lig_matches = [res for res in modeller.topology.residues() if res.name.upper() == keep_resname]
    if len(lig_matches) == 0:
        raise RuntimeError(
            f"No ligand residue with resname '{keep_resname}' found. "
            f"Please check HETATM lines in the input PDB."
        )

    # 选择要保留的 ligand（默认只保留一个，避免 FEP 误 alchemize 两个）
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
    # 3) 构建删除列表：删水、离子、以及所有“非标准蛋白/核酸/保留ligand”
    # -----------------------
    delete_res = []
    for res in modeller.topology.residues():
        name = res.name.upper()

        # 保留目标 ligand（指定那一个/些）
        if name == keep_resname and res in keep_lig_set:
            continue

        # 删同 resname 的其它 ligand 副本
        if name == keep_resname and res not in keep_lig_set:
            delete_res.append(res)
            continue

        # 水
        if name in WATER_NAMES:
            if not args.keep_water:
                delete_res.append(res)
            continue

        # 离子
        if name in COMMON_IONS:
            if not args.keep_ions:
                delete_res.append(res)
            continue

        # 标准蛋白残基：保留
        if name in STANDARD_AA:
            continue

        # 核酸：按需保留
        if args.keep_nucleic_acids and name in STANDARD_NA:
            continue

        # 其它全删（这一步就是“真正干净”的关键）
        delete_res.append(res)

    if delete_res:
        modeller.delete(delete_res)

    # -----------------------
    # 4) 写出
    # -----------------------
    with open(args.out_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)

    # -----------------------
    # 5) 打印 summary（你可以用它确认删干净了）
    # -----------------------
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
