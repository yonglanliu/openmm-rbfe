## protein_process.py
This utility cleans and standardizes a protein–ligand PDB for downstream MD/FEP workflows:
* Fixes common PDB issues with PDBFixer (nonstandard residues, missing atoms, hydrogens).
* Keeps only the specified ligand by residue name (--keep_resname).
* Removes all other heterogens robustly (waters/ions optional).
* Optionally restricts ligand selection by chain and/or residue id to avoid keeping the wrong copy when multiple ligands exist.

**Keep ligand by resname:**
```bash
python tools/fix_pdb_keep_ligand.py \
        -i data/pdb/[PDBid].pdb \
        -o data/protein/[PDBid]_fixed_keep[ligand Resname].pdb \
        --keep_resname [ligand Resname] \
```  
---
**Choose the correct ligand when multiple copies exist:**
```bash
python tools/fix_pdb_keep_ligand.py \
        -i data/pdb/[PDBid].pdb \
        -o data/protein/[PDBid]_fixed_keep[ligand Resname].pdb \
        --keep_resname [ligand Resname] \
        --keep_chain [chain id] \
        --keep_resid [residue id]
```
---
**Arguments:**
* <code>-i, --in_pdb</code> (required): Input PDB file.
* <code>-o, --out_pdb</code>  (required): Output cleaned PDB file.
* <code>--keep_resname</code> (required): Ligand residue name to keep (e.g., LIG, BNZ).
* <code>--keep_chain</code>  (optional): Keep ligand only on this chain (e.g., A).
* <code>--keep_resid</code>  (optional): Keep ligand only with this residue id (e.g., 401).
* <code>--ph</code>  (optional, default 7.0): pH used when adding hydrogens.
* <code>--keep_water</code> : Keep crystallographic waters (default: remove).
* <code>--keep_ions</code> : Keep common ions (default: remove).
* <code>--keep_nucleic_acids</code> : Keep DNA/RNA residues (default: remove).
* <code>--keep_all_matching_ligands</code> : Keep all ligands matching --keep_resname (default: keep only one).
---

**Output and sanity checks**
* After running, the script prints a summary like:
* Ligands found before
* Ligands kept after
* “Other heterogens after” (should be empty for a clean protein–ligand system)
```bash
Input : input.pdb
Output: cleaned.pdb
Keep ligand resname: LIG
Ligands found BEFORE: [('LIG','A','401'), ('LIG','B','402')]
Ligands kept AFTER : [('LIG','A','401')]
Other heterogens AFTER (should be empty): []
```
---
**Notes / Best practices:**
* If you see No ligand residue with resname 'XXX' found, check the input PDB HETATM residue name (columns 18–20) and use that exact resname in --keep_resname.
* When multiple ligands exist, strongly recommend specifying --keep_chain and/or --keep_resid to avoid keeping the wrong ligand.
* This script intentionally removes all non-protein/non-selected-ligand heterogens by default to produce a “clean” PDB suitable for automated parameterization and system building.

### extract_ligand_to_sdf.py
This utility extracts the ligand 3D coordinates from a protein–ligand complex PDB and writes an SDF with correct chemistry (bond orders, valence, explicit H) using RDKit.
**Typical use case:**
* You have a cleaned complex PDB (e.g., from fix_pdb_keep_ligand.py),
* You want a ligand SDF that preserves the pose from the PDB while restoring chemical connectivity from a trusted SMILES.
---
**Basic command:**
Extract ligand coordinates by residue name and write an SDF:
```bash
python tools/extract_ligand_pose_to_sdf.py \
        -i data/protein/fixed_complex.pdb \
        --resname [ligand_resname] \
        --smiles "c1ccccc1" \
        -o data/ligands/[ligand_resnamd].sdf
```
---

**Arguments:**
* <code>-i, --pdb</code> (required): Input complex PDB (protein + ligand).(Recommended: use a cleaned PDB where only the target ligand remains.)
* <code>--resname</code>  (required): Ligand residue name in the PDB (e.g., BNZ, MBN, LIG).
* <code>--smiles</code>  (required): Ligand SMILES (used to build correct bond orders/chemistry).
* <code>-o, --out_sdf</code>  (required): Output SDF path.
* <code>--out_lig_pdb</code>  (optional): Also write a ligand-only PDB exported from the RDKit molecule.
---

**Notes / Best practices:**
* Heavy atom count must match between:
  * ligand atoms extracted from the PDB (by --resname), and
  * the heavy atoms in the SMILES. If not, the script will error with a “Heavy atom count mismatch”.
* The atom mapping uses “attached hydrogen count” grouping; it works very well for **simple aromatics (e.g., benzene/toluene)** and other ligands with clear H-count patterns.
* For more complex ligands (symmetry, ambiguous H placement, missing H coordinates), mapping may fail.
* If your PDB does not contain ligand hydrogens or the mapping fails:
  * consider providing an SDF from RCSB/ligand expo as a reference, or
  * regenerate a ligand PDB/SDF with consistent atom ordering (then align/transfer coordinates).
---