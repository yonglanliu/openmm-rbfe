import os
import sys

import pandas as pd
import streamlit as st

# -----------------------------
# Path setup so we can import tools/
# -----------------------------
VIZ_DIR = os.path.dirname(os.path.abspath(__file__))      # .../openmm-rbfe/viz
ROOT = os.path.dirname(VIZ_DIR)                           # .../openmm-rbfe
sys.path.insert(0, ROOT)

from tools.summarize_from_qc import summarize_ligand

# -----------------------------
# Streamlit page config
# -----------------------------
st.set_page_config(
    page_title="RBFE QC Summary",
    page_icon="📊",
    layout="wide",
)

# -----------------------------
# Session state: store summaries
# -----------------------------
# Each entry will be keyed by (lig_tag, leg)
# and hold: rows, mean, std, sem, out_csv
if "summaries" not in st.session_state:
    st.session_state["summaries"] = {}  # {(lig_tag, leg): {...}}


# -----------------------------
# Sidebar controls
# -----------------------------
st.sidebar.title("RBFE QC Summary")

st.sidebar.markdown(
    """
Choose the **leg** and **ligand**, then click
**Summarize QC Results** to aggregate MBAR outputs across seeds.
You can run **complex** and **hydration** one after another;
both results will stay on this page.
"""
)

leg = st.sidebar.selectbox("Leg", ["complex", "hydration"])
lig_tag = st.sidebar.text_input("Ligand tag (e.g. MBN, BNZ)", "")

run_button = st.sidebar.button("Summarize QC Results", type="primary")

st.sidebar.markdown("---")
st.sidebar.markdown(
    """
**Expected QC files**

- `results/qc/<LIG>_complex_seedXXXX_mbar_qc.csv`  
- `results/qc/<LIG>_hydration_seedXXXX_mbar_qc.csv`

Summary CSVs are written to `results/`.
"""
)

# -----------------------------
# Main content
# -----------------------------
st.title("📊 RBFE QC Summary Dashboard")

st.markdown(
    """
This tool summarizes **MBAR QC results** across multiple seeds
for a given ligand and leg (**complex** or **hydration**).

You can compute multiple summaries (e.g., complex and hydration for the same ligand)
and they will all remain visible below.
"""
)

# -----------------------------
# Handle button: compute and store in session_state
# -----------------------------
if run_button:
    if not lig_tag:
        st.error("Please enter a ligand tag in the sidebar.")
    else:
        try:
            rows, mean, std, sem, out_csv = summarize_ligand(lig_tag, leg=leg)

            key = (lig_tag, leg)
            st.session_state["summaries"][key] = {
                "rows": rows,
                "mean": mean,
                "std": std,
                "sem": sem,
                "out_csv": out_csv,
            }

            st.success(f"{leg.capitalize()} summary for {lig_tag} created and stored.")

        except FileNotFoundError as e:
            st.error(str(e))
        except Exception as e:
            st.error(f"Unexpected error: {e}")

# -----------------------------
# Render all stored summaries
# -----------------------------
if st.session_state["summaries"]:
    st.markdown("## Stored summaries")

    # One expandable block per (ligand, leg)
    for (lig, leg_name), data in st.session_state["summaries"].items():
        rows = data["rows"]
        mean = data["mean"]
        std = data["std"]
        sem = data["sem"]
        out_csv = data["out_csv"]

        with st.expander(f"{lig} — {leg_name} leg", expanded=True):
            # Top metrics
            col1, col2, col3 = st.columns(3)
            col1.metric("Mean ΔG (kJ/mol)", f"{mean:.3f}")
            col2.metric("Std (kJ/mol)", f"{std:.3f}")
            col3.metric("SEM (kJ/mol)", f"{sem:.3f}")

            # Per-seed DataFrame
            df = pd.DataFrame(
                {
                    "seed": [r[0] for r in rows],
                    "dG_kJmol": [r[1] for r in rows],
                    "SEM_kJmol": [r[2] for r in rows],
                }
            ).set_index("seed")

            st.markdown("**Per-seed results**")
            st.dataframe(
                df.style.format({"dG_kJmol": "{:.3f}", "SEM_kJmol": "{:.3f}"}),
                use_container_width=True,
            )

            # Small bar chart
            if len(df) > 1:
                st.markdown("**ΔG per seed**")
                st.bar_chart(df["dG_kJmol"])

            # Download button
            st.markdown("**Export**")
            st.caption(f"Summary CSV: `{out_csv}`")
            try:
                with open(out_csv, "rb") as f:
                    st.download_button(
                        label=f"Download {leg_name} summary CSV for {lig}",
                        data=f,
                        file_name=os.path.basename(out_csv),
                        mime="text/csv",
                        key=f"download_{lig}_{leg_name}",
                    )
            except FileNotFoundError:
                st.warning("Summary CSV file not found on disk (maybe moved or deleted).")
else:
    st.info(
        "Use the controls in the sidebar to select leg and ligand, "
        "then click **Summarize QC Results**. All computed summaries will appear here."
    )
