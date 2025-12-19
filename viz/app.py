import sys
from pathlib import Path

import pandas as pd
import streamlit as st

# -------------------------------------------------------------------
# Path setup so we can import tools/ from the project root
# -------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]   # .../openmm-rbfe
sys.path.insert(0, str(ROOT))

from tools.summarize_from_qc import summarize_ligand  # must accept leg="complex"/"hydration"

# -------------------------------------------------------------------
# Streamlit page config
# -------------------------------------------------------------------
st.set_page_config(page_title="OpenMM FEP Dashboard", layout="wide")

# -------------------------------------------------------------------
# Original header content (as you had it)
# -------------------------------------------------------------------
st.title("OpenMM FEP Dashboard")
st.caption("Hydration FEP (solvent leg) — results, repeats, and λ-space QC")

root = Path(".").resolve()
results_dir = root / "results"
qc_dir = results_dir / "qc"

if not results_dir.exists():
    st.error("No ./results directory found. Run the pipeline first.")
    st.stop()

st.markdown("### What this dashboard reads")
st.markdown(
    "- `results/ddg_summary_repeats.csv`\n"
    "- `results/qc/*_mbar_qc.csv` (per ligand & seed)\n"
)

st.success("Use the left sidebar to navigate pages: Overview / QC Explorer / Window Inspector")

# -------------------------------------------------------------------
# Session state: store multiple summaries (complex / hydration, ligands)
# -------------------------------------------------------------------
# Key: (lig_tag, leg) -> dict(rows, mean, std, sem, out_csv)
if "summaries" not in st.session_state:
    st.session_state["summaries"] = {}

# -------------------------------------------------------------------
# Sidebar controls for summarization
# -------------------------------------------------------------------
st.sidebar.title("RBFE QC Summary")

st.sidebar.markdown(
    """
Choose the **leg** and **ligand**, then click
**Summarize QC Results** to aggregate MBAR outputs across seeds.

You can run **complex** and **hydration** for the same ligand
(one after another) and see both below on this page.
"""
)

leg = st.sidebar.selectbox("Leg", ["complex", "hydration"])
lig_tag = st.sidebar.text_input("Ligand tag (e.g. MBN, BNZ)", "")

run_button = st.sidebar.button("Summarize QC Results", type="primary")

st.sidebar.markdown("---")
st.sidebar.markdown(
    f"""
**Expected QC files**

- `{qc_dir}/<LIG>_complex_seedXXXX_mbar_qc.csv`  
- `{qc_dir}/<LIG>_hydration_seedXXXX_mbar_qc.csv`

Summary CSVs are written to `{results_dir}/`.
"""
)

# -------------------------------------------------------------------
# Interactive summary section
# -------------------------------------------------------------------
st.markdown("## Interactive QC Summary")

st.markdown(
    """
This section summarizes **MBAR QC results** across independent seeds
for a given ligand and leg (**complex** or **hydration**).

Each (ligand, leg) summary you compute is stored, so you can
compare complex vs hydration or multiple ligands on the same page.
"""
)

# Handle button: compute summary and store in session_state
if run_button:
    if not lig_tag:
        st.error("Please enter a ligand tag in the sidebar.")
    else:
        try:
            rows, mean, std, sem, out_csv = summarize_ligand(
                lig_tag,
                qc_dir=str(qc_dir),
                outdir=str(results_dir),
                leg=leg,
            )

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

# -------------------------------------------------------------------
# Render all stored summaries
# -------------------------------------------------------------------
if st.session_state["summaries"]:
    st.markdown("### Stored summaries")

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

            # Bar chart: ΔG per seed
            if len(df) > 1:
                st.markdown("**ΔG per seed**")
                st.bar_chart(df["dG_kJmol"])

            # Cumulative ΔG vs window index, per seed
            st.markdown("**Cumulative ΔG across λ windows (per seed)**")

            all_cum = []
            for seed in df.index:
                qc_path = qc_dir / f"{lig}_{leg_name}_seed{seed}_mbar_qc.csv"
                if not qc_path.exists():
                    st.warning(f"QC file not found for seed {seed}: {qc_path}")
                    continue

                try:
                    qc_df = pd.read_csv(qc_path)
                    if "i" not in qc_df.columns or "dG_cum_to_i_kJmol" not in qc_df.columns:
                        st.warning(f"Unexpected QC format for seed {seed}: {qc_path}")
                        continue

                    df_seed = qc_df[["i", "dG_cum_to_i_kJmol"]].copy()
                    df_seed["seed"] = seed
                    all_cum.append(df_seed)
                except Exception as e:
                    st.warning(f"Failed to read QC for seed {seed}: {e}")

            if all_cum:
                cum_df = pd.concat(all_cum, ignore_index=True)
                plot_df = cum_df.pivot_table(
                    index="i",
                    columns="seed",
                    values="dG_cum_to_i_kJmol",
                ).sort_index()

                st.line_chart(plot_df)
            else:
                st.info("No cumulative QC data available to plot for this summary.")

            # Download button
            st.markdown("**Export**")
            st.caption(f"Summary CSV: `{out_csv}`")
            try:
                with open(out_csv, "rb") as f:
                    st.download_button(
                        label=f"Download {leg_name} summary CSV for {lig}",
                        data=f,
                        file_name=Path(out_csv).name,
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
