import streamlit as st
from pathlib import Path
from viz.io import scan_qc_csv, available_ligands, available_seeds, qc_path_for, load_qc_csv
from viz.plots import find_worst_windows

st.title("Window Inspector")

results_dir = Path("results")
qc_items = scan_qc_csv(results_dir)

if not qc_items:
    st.error("No QC CSVs found in results/qc. Generate QC first.")
    st.stop()

ligands = available_ligands(qc_items)
ligand = st.sidebar.selectbox("Ligand", ligands)
seeds = available_seeds(qc_items, ligand)
seed = st.sidebar.selectbox("Seed", seeds)

p = qc_path_for(qc_items, ligand, seed)
df = load_qc_csv(p)

st.caption(f"Loaded: {p}")

st.subheader("Worst windows (by |ΔG_adj|)")
worst_mag, worst_err = find_worst_windows(df, topn=5)
st.dataframe(worst_mag[["i","le_i","ls_i","le_ip1","ls_ip1","dG_adj_kJmol","dG_adj_err_kJmol"]], use_container_width=True)

st.subheader("Worst windows (by uncertainty)")
st.dataframe(worst_err[["i","le_i","ls_i","le_ip1","ls_ip1","dG_adj_kJmol","dG_adj_err_kJmol"]], use_container_width=True)

st.info(
    "Interpretation:\n"
    "- Large |ΔG_adj| or large uncertainty usually appears near sterics decoupling close to λ_sterics→0.\n"
    "- If spikes occur only for one seed, it suggests sampling noise (increase prod steps for that region).\n"
    "- If spikes occur for all seeds in the same region, densify windows there (add more λ points)."
)
