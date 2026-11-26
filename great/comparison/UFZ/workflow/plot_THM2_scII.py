"""Plots of the UFZ Scenario II comparison files (THM2, Load A), read from the CSV files written by
export_THM2_scII.py (GREAT Cell Modelling Specification v2, Section 8.2).

Figures (written to <results_dir>):
  UFZ_THM2_scII_histories.png : T, p, eps_vol at P0...P270 and inlet p_in, T_in, Q, V_cum
  UFZ_THM2_scII_profiles.png  : T(theta), p(theta), eps_vol(theta) on r_p = 0.065 m at the snapshot times

Optional --ref <dir> adds the four-point histories of OGS model v2 (with rubber sheath,
files <case>_points.csv of ../../THM_model_v2/results) at P0 for comparison.

Usage: python plot_THM2_scII.py <results_dir> [--ref ../../THM_model_v2/results]
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("res")
ap.add_argument("--cases", nargs="+", default=["THM2-GW-A-T300", "THM2-GN-A-T300"])
ap.add_argument("--ref", default=None, help="folder with model v2 results (with sheath)")
a = ap.parse_args()

NAME = {"GW": "Greywacke", "GN": "Gneiss"}
PCOL = {"P0": "tab:blue", "P90": "tab:orange", "P180": "tab:green", "P270": "tab:red"}
SNAP = [2000, 2060, 3800, 5600]


def shade(ax):
    ax.axvspan(2000, 2120, color="0.9", zorder=0)
    ax.axvline(1800, color="0.6", lw=0.6, ls=":")
    ax.set_xlim(0, 5600); ax.set_xlabel("t [s]"); ax.grid(alpha=0.3)


# ------------------------------------------------------------------ histories
fig, ax = plt.subplots(len(a.cases), 4, figsize=(20, 4.6 * len(a.cases)), squeeze=False)
for row, c in enumerate(a.cases):
    rock = c.split("-")[1]
    P = pd.read_csv(f"{a.res}/UFZ_{c}_points.csv")
    I = pd.read_csv(f"{a.res}/UFZ_{c}_inlet.csv")
    for n, col in PCOL.items():
        ls = "-" if n in ("P0", "P90") else (0, (4, 3))      # P180 = P0, P270 = P90 by symmetry
        ax[row, 0].plot(P.t, P[f"T_{n}"], color=col, ls=ls, label=n)
        ax[row, 1].plot(P.t, P[f"p_{n}"], color=col, ls=ls, label=n)
        ax[row, 2].plot(P.t, P[f"eps_vol_{n}"] * 1e6, color=col, ls=ls, label=n)
    if a.ref and os.path.exists(f"{a.ref}/{c}_points.csv"):
        R = pd.read_csv(f"{a.ref}/{c}_points.csv")
        kw = dict(color="k", lw=0.9, ls=":", label="P0, model v2 (with sheath)")
        ax[row, 0].plot(R.t_s, R.T_P0_degC, **kw)
        ax[row, 1].plot(R.t_s, R.p_P0_MPa, **kw)
        ax[row, 2].plot(R.t_s, R.epsvol_P0 * 1e6, **kw)
    a4 = ax[row, 3]
    a4.plot(I.t, I.p_in, "k-", label=r"$p_\mathrm{in}$ [MPa]")
    a4.plot(I.t, I.T_in / 10, "r-", label=r"$T_\mathrm{in}$ / 10 [°C]")
    a4.set_ylabel(r"$p_\mathrm{in}$ [MPa], $T_\mathrm{in}/10$ [°C]")
    a5 = a4.twinx()
    a5.plot(I.t, I.Q_ml_min * 10, "g--", label="Q × 10 [ml/min]")
    a5.plot(I.t, I.V_cum_ml, "b--", label=r"$V_\mathrm{cum}$ [ml]")
    a5.set_ylabel("Q × 10 [ml/min], V [ml]"); a5.set_ylim(0, 7)
    h1, l1 = a4.get_legend_handles_labels(); h2, l2 = a5.get_legend_handles_labels()
    a4.legend(h1 + h2, l1 + l2, fontsize=7, loc="center right")
    ax[row, 0].set_ylabel(f"{NAME[rock]}\nT [°C]")
    ax[row, 1].set_ylabel("p (gauge) [MPa]")
    ax[row, 2].set_ylabel(r"$\varepsilon_\mathrm{vol}$ [$\mu\varepsilon$]")
    for x in ax[row]:
        shade(x)
    for x in ax[row, :3]:
        x.legend(fontsize=7)
    ax[row, 0].set_title(f"{c}: temperature on r = 0.065 m")
    ax[row, 1].set_title(f"{c}: pore pressure on r = 0.065 m")
    ax[row, 2].set_title(f"{c}: volumetric strain on r = 0.065 m")
    ax[row, 3].set_title(f"{c}: inlet (centre node)")
fig.suptitle("UFZ / OGS 6.5.9 - GREAT cell Scenario II, THM2, Load A (no rubber sheath)", fontsize=13)
fig.tight_layout()
fig.savefig(f"{a.res}/UFZ_THM2_scII_histories.png", dpi=150)

# ------------------------------------------------------------------ circumferential profiles
fig, ax = plt.subplots(len(a.cases), 3, figsize=(16, 4.6 * len(a.cases)), squeeze=False)
cmap = plt.get_cmap("viridis")
for row, c in enumerate(a.cases):
    rock = c.split("-")[1]
    D = pd.read_csv(f"{a.res}/UFZ_{c}_profiles.csv")
    for i, ts in enumerate(SNAP):
        d = D[np.isclose(D.t, ts)]
        kw = dict(color=cmap(i / (len(SNAP) - 1)), label=f"t = {ts} s")
        ax[row, 0].plot(d.theta_deg, d["T"], **kw)
        ax[row, 1].plot(d.theta_deg, d["p"], **kw)
        ax[row, 2].plot(d.theta_deg, d.eps_vol * 1e6, **kw)
    ax[row, 0].set_ylabel(f"{NAME[rock]}\nT [°C]")
    ax[row, 1].set_ylabel("p (gauge) [MPa]")
    ax[row, 2].set_ylabel(r"$\varepsilon_\mathrm{vol}$ [$\mu\varepsilon$]")
    for x in ax[row]:
        x.set_xlim(0, 360); x.set_xticks(range(0, 361, 45)); x.set_xlabel(r"$\theta$ [deg]")
        x.grid(alpha=0.3); x.legend(fontsize=7)
    ax[row, 0].set_title(f"{c}: T on r = 0.065 m")
    ax[row, 1].set_title(f"{c}: p on r = 0.065 m")
    ax[row, 2].set_title(rf"{c}: $\varepsilon_{{vol}}$ on r = 0.065 m")
fig.tight_layout()
fig.savefig(f"{a.res}/UFZ_THM2_scII_profiles.png", dpi=150)
print(f"saved {a.res}/UFZ_THM2_scII_histories.png, UFZ_THM2_scII_profiles.png")
