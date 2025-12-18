import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import streamlit as st
from pathlib import Path

st.set_page_config(page_title="OpenMM FEP Dashboard", layout="wide")

st.title("OpenMM FEP Dashboard")
st.caption("Hydration FEP (solvent leg) — results, repeats, and λ-space QC")

root = Path(".").resolve()
results_dir = root / "results"

if not results_dir.exists():
    st.error("No ./results directory found. Run the pipeline first.")
    st.stop()

st.markdown("### What this dashboard reads")
st.markdown(
    "- `results/ddg_summary_repeats.csv`\n"
    "- `results/qc/*_mbar_qc.csv` (per ligand & seed)\n"
)

st.success("Use the left sidebar to navigate pages: Overview / QC Explorer / Window Inspector")
