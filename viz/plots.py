from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_adjacent(df: pd.DataFrame, title: str):
    x = df["i"].to_numpy()
    y = df["dG_adj_kJmol"].to_numpy()
    yerr = df["dG_adj_err_kJmol"].to_numpy()

    fig = plt.figure()
    plt.errorbar(x[:-1], y[:-1], yerr=yerr[:-1], fmt="o")  # last row may be blank for adjacent
    plt.xlabel("Adjacent window index i (i → i+1)")
    plt.ylabel("ΔG_adj (kJ/mol)")
    plt.title(title)
    plt.tight_layout()
    return fig

def plot_cumulative(df: pd.DataFrame, title: str):
    x = df["i"].to_numpy()
    y = df["dG_cum_to_i_kJmol"].to_numpy()
    yerr = df["dG_cum_err_to_i_kJmol"].to_numpy()

    fig = plt.figure()
    plt.errorbar(x, y, yerr=yerr, fmt="o")
    plt.xlabel("Window index i (state 0 → i)")
    plt.ylabel("Cumulative ΔG (kJ/mol)")
    plt.title(title)
    plt.tight_layout()
    return fig

def plot_overlay_adjacent(dfs: list[tuple[int, pd.DataFrame]], title: str):
    fig = plt.figure()
    for seed, df in dfs:
        x = df["i"].to_numpy()
        y = df["dG_adj_kJmol"].to_numpy()
        plt.plot(x[:-1], y[:-1], marker="o", alpha=0.7, label=f"seed {seed}")
    plt.xlabel("Adjacent window index i (i → i+1)")
    plt.ylabel("ΔG_adj (kJ/mol)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    return fig

def find_worst_windows(df: pd.DataFrame, topn: int = 5):
    # rank by absolute adjacent contribution and by uncertainty
    tmp = df.copy()
    tmp = tmp.iloc[:-1].copy()  # last row may not have adjacent
    tmp["abs_dG_adj"] = np.abs(tmp["dG_adj_kJmol"])
    tmp["abs_err"] = np.abs(tmp["dG_adj_err_kJmol"])

    worst_by_mag = tmp.sort_values("abs_dG_adj", ascending=False).head(topn)
    worst_by_err = tmp.sort_values("abs_err", ascending=False).head(topn)
    return worst_by_mag, worst_by_err

import plotly.graph_objects as go

def plot_adjacent_plotly(df, title: str):
    x = df["i"][:-1]
    y = df["dG_adj_kJmol"][:-1]
    yerr = df["dG_adj_err_kJmol"][:-1]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers+lines",
            error_y=dict(type="data", array=yerr, visible=True),
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title="Adjacent window index i (i → i+1)",
        yaxis_title="ΔG_adj (kJ/mol)",
        hovermode="closest",
    )
    return fig


def plot_cumulative_plotly(df, title: str):
    x = df["i"]
    y = df["dG_cum_to_i_kJmol"]
    yerr = df["dG_cum_err_to_i_kJmol"]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers+lines",
            error_y=dict(type="data", array=yerr, visible=True),
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title="Window index i (state 0 → i)",
        yaxis_title="Cumulative ΔG (kJ/mol)",
        hovermode="closest",
    )
    return fig

def plot_adjacent_plotly_highlight(df, title: str, topn: int = 5):
    """
    Adjacent ΔG plot with automatic highlighting of worst windows:
      - largest |ΔG_adj|
      - largest uncertainty
    """
    # adjacent data (last row may not have adjacent)
    dfa = df.iloc[:-1].copy()

    x = dfa["i"].to_numpy()
    y = dfa["dG_adj_kJmol"].to_numpy()
    yerr = dfa["dG_adj_err_kJmol"].to_numpy()

    # indices for highlights
    idx_mag = np.argsort(np.abs(y))[::-1][:topn]
    idx_err = np.argsort(np.abs(yerr))[::-1][:topn]

    fig = go.Figure()

    # base trace
    fig.add_trace(
        go.Scatter(
            x=x, y=y,
            mode="markers+lines",
            error_y=dict(type="data", array=yerr, visible=True),
            name="ΔG_adj",
        )
    )

    # highlight: largest magnitude
    fig.add_trace(
        go.Scatter(
            x=x[idx_mag], y=y[idx_mag],
            mode="markers",
            marker=dict(symbol="circle", size=12),
            name=f"Top {topn} |ΔG_adj|",
        )
    )

    # highlight: largest uncertainty
    fig.add_trace(
        go.Scatter(
            x=x[idx_err], y=y[idx_err],
            mode="markers",
            marker=dict(symbol="x", size=12),
            name=f"Top {topn} uncertainty",
        )
    )

    # annotate with (λ_elec, λ_sterics) for context (optional but useful)
    for ii in np.unique(np.concatenate([idx_mag, idx_err])):
        le_i = float(dfa.loc[dfa.index[ii], "le_i"])
        ls_i = float(dfa.loc[dfa.index[ii], "ls_i"])
        fig.add_annotation(
            x=x[ii], y=y[ii],
            text=f"i={int(x[ii])}<br>le={le_i:.2f}, ls={ls_i:.2f}",
            showarrow=True, arrowhead=2, ax=20, ay=-30,
        )

    fig.update_layout(
        title=title,
        xaxis_title="Adjacent window index i (i → i+1)",
        yaxis_title="ΔG_adj (kJ/mol)",
        hovermode="closest",
    )
    return fig