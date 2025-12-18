import streamlit as st
from pathlib import Path
from viz.io import find_ddg_summary, load_ddg_summary, scan_qc_csv, available_ligands

st.title("Overview")

results_dir = Path("results")
ddg_path = find_ddg_summary(results_dir)

if ddg_path is None:
    st.warning("Missing results/ddg_summary_repeats.csv (run repeats first).")
else:
    df = load_ddg_summary(ddg_path)
    st.subheader("ΔΔG repeats summary (CSV)")
    st.dataframe(df, use_container_width=True)

qc_items = scan_qc_csv(results_dir)
st.subheader("QC availability")
if not qc_items:
    st.warning("No QC CSVs found in results/qc. Generate QC via save_mbar_qc().")
else:
    ligs = available_ligands(qc_items)
    st.write(f"Found QC for ligands: {', '.join(ligs)}")
    st.write(f"Total QC files: {len(qc_items)}")
