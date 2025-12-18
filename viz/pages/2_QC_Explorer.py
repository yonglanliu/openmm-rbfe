import streamlit as st
from pathlib import Path
from viz.io import scan_qc_csv, available_ligands, available_seeds, qc_path_for, load_qc_csv
from viz.plots import plot_adjacent_plotly_highlight, plot_cumulative_plotly, plot_overlay_adjacent


st.title("QC Explorer (λ-space)")

results_dir = Path("results")
qc_items = scan_qc_csv(results_dir)

if not qc_items:
    st.error("No QC CSVs found in results/qc. Generate QC first.")
    st.stop()

ligands = available_ligands(qc_items)
ligand = st.sidebar.selectbox("Ligand", ligands)

seeds = available_seeds(qc_items, ligand)
mode = st.sidebar.radio("Mode", ["Single seed", "Overlay seeds"])

if mode == "Single seed":
    seed = st.sidebar.selectbox("Seed", seeds)
    p = qc_path_for(qc_items, ligand, seed)
    df = load_qc_csv(p)

    st.caption(f"Loaded: {p}")
    

    st.plotly_chart(
    plot_adjacent_plotly_highlight(df, f"{ligand} — Adjacent ΔG (seed {seed})", topn=5),
    use_container_width=True,
)


    st.plotly_chart(
        plot_cumulative_plotly(df, f"{ligand} — Cumulative ΔG (seed {seed})"),
        use_container_width=True,
    )

    st.dataframe(df, use_container_width=True)

else:
    selected = st.sidebar.multiselect("Seeds to overlay", seeds, default=seeds[: min(3, len(seeds))])
    dfs = []
    for s in selected:
        p = qc_path_for(qc_items, ligand, s)
        if p is None:
            continue
        dfs.append((s, load_qc_csv(p)))

    if len(dfs) < 1:
        st.warning("Select at least one seed.")
    else:
        st.plotly_chart(plot_overlay_adjacent(dfs, f"{ligand} — Adjacent ΔG overlay"))
        st.info("Tip: If one seed shows a spike not seen in others, that window likely needs more sampling or denser spacing.")
