"""Overview plot of the GREAT cell THM model (Load A), with reference curves from the
TH model (centre temperature) and the M1 model (volumetric strain on the probe circle).

Usage: python plot_THM_overview.py [out_dir] [png] [--th TH_out_dir] [--m1 M1_out_dir]
"""
import argparse
import glob
import re
import os
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("out", nargs="?", default="out")
ap.add_argument("png", nargs="?", default="THM_loadA_overview.png")
ap.add_argument("--load", default="A")
ap.add_argument("--th", default=None, help="TH model output folder (centre temperature reference)")
ap.add_argument("--m1", default=None, help="M1 model output folder (strain reference)")
a = ap.parse_args()

P0, R, K, MU, Q = 1e6, 0.1, 1.409e-14, 1e-3, 1.6666667e-5
ROCKS = {"greywacke": "tab:green", "gneiss": "tab:red"}


def series(prefix):
    fs = glob.glob(f"{prefix}_ts_*_t_*.vtu")
    t = np.array([float(re.search(r"_t_([0-9.]+)\.vtu$", f).group(1)) for f in fs])
    o = np.argsort(t)
    return t[o], [fs[i] for i in o]


def corner_nodes(m):
    """T and p are linear (order 1): values are only defined on the corner nodes of the
    quadratic triangles; mid-edge nodes carry zeros in the output and are skipped"""
    c = m.cells_dict[22] if 22 in m.cells_dict else m.cells_dict[5]
    return np.unique(c[:, :3])


def probe(f):
    m = pv.read(f)
    th = np.degrees(np.arctan2(m.points[:, 1], m.points[:, 0])) % 360
    e = m.point_data["epsilon"]
    o = np.argsort(th)
    return th[o], (e[o, 0] + e[o, 1]) * 1e4      # [1e-2 %]


fig, ax = plt.subplots(2, 3, figsize=(17, 10))
(a1, a2, a3), (b1, b2, b3) = ax

# --- centre node: T and p over time --------------------------------------------------------
for rock, c in ROCKS.items():
    t, fs = series(f"{a.out}/THM_{rock}_load{a.load}_center")
    T = np.array([pv.read(f).point_data["temperature"][0] for f in fs])
    p = np.array([pv.read(f).point_data["pressure"][0] for f in fs])
    a1.plot(t, T, color=c, label=f"T, THM {rock}")
    a2.plot(t, (p - P0) / 1e6, color=c, lw=2.5 if rock == "greywacke" else 1.2,
            label=f"p - p0, THM {rock}")
    if rock == "greywacke":
        a2b = a2.twinx()
    a2b.plot(t, T, color=c, ls=":", lw=1)
if a.th:
    t, fs = series(f"{a.th}/TH_granite_center")
    T = np.array([pv.read(f).point_data["temperature"][0] for f in fs])
    a1.plot(t, T, "k--", lw=1, label="T, TH model (reference)")
a1.set_xlabel("t [s]"); a1.set_ylabel("T [K]"); a1.set_title("Centre node: temperature")
a1.axvspan(4000, 4050, color="0.85", zorder=0); a1.legend(fontsize=8, loc="lower right")
a2.set_xlim(3990, 4100); a2.set_xlabel("t [s]"); a2.set_ylabel("p - p0 [MPa] (solid)")
a2b.set_ylabel("T [K] (dotted)")
a2.axvspan(4000, 4050, color="0.85", zorder=0)
a2.set_title("Centre node, zoom on injection"); a2.legend(fontsize=8, loc="center right")

# --- radial profiles (greywacke) ---------------------------------------------------------------
rock = "greywacke"
t_dom, f_dom = series(f"{a.out}/THM_{rock}_load{a.load}_domain")
dom = {round(t, 3): f for t, f in zip(t_dom, f_dom)}
cmap = plt.get_cmap("viridis")
times_T = [4000, 4010, 4025, 4050, 4100, 4200, 4500, 5000]
for i, t in enumerate(times_T):
    m = pv.read(dom[t])
    cn = corner_nodes(m)
    r = np.hypot(m.points[cn, 0], m.points[cn, 1])
    a3.plot(r, m.point_data["temperature"][cn], ".", ms=1.5, color=cmap(i / (len(times_T) - 1)),
            label=f"t = {t} s")
a3.set_xlabel("r [m]"); a3.set_ylabel("T [K]"); a3.set_title(f"Radial temperature profiles ({rock})")
a3.legend(markerscale=6, fontsize=8, loc="lower right")

times_p = [4000, 4001, 4025, 4050, 4060]
for i, t in enumerate(times_p):
    m = pv.read(dom[t])
    cn = corner_nodes(m)
    r = np.hypot(m.points[cn, 0], m.points[cn, 1])
    s = r > 0
    b1.plot(r[s], (m.point_data["pressure"][cn][s] - P0) / 1e6, "o", ms=6 - i, mfc="none",
            color=cmap(i / (len(times_p) - 1)), ls="none", label=f"t = {t} s")
rr = np.logspace(-3, np.log10(R), 200)
b1.plot(rr, Q * MU / (2 * np.pi * K) * np.log(R / rr) / 1e6, "k-", lw=0.8, label="Thiem (steady)")
b1.set_xscale("log"); b1.set_xlabel("r [m]"); b1.set_ylabel("p - p0 [MPa]")
b1.set_title(f"Radial pressure profiles ({rock})"); b1.legend(fontsize=8)

# --- volumetric strain on the probe circle -----------------------------------------------------
times_e = [5, 4000, 4050, 5000]
ls = dict(zip(times_e, ["-", "--", "-.", ":"]))
for rock, c in ROCKS.items():
    t, fs = series(f"{a.out}/THM_{rock}_load{a.load}_probe_circle")
    pf = {round(x, 3): f for x, f in zip(t, fs)}
    for te in times_e:
        th, ev = probe(pf[te])
        b2.plot(th, ev, ls[te], color=c, lw=1.2, label=f"{rock}, t = {te} s")
    if a.m1:
        th, ev = probe(f"{a.m1}/M1_{rock}_load{a.load}_probe_circle_ts_1_t_1.000000.vtu")
        b2.plot(th[::4], ev[::4], "o", ms=3, mfc="none", color="k",
                label="M1 model" if rock == "greywacke" else None)
    # mean volumetric strain on the probe circle vs time
    tt, em = [], []
    for x, f in zip(t, fs):
        e = pv.read(f).point_data["epsilon"]
        tt.append(x); em.append(np.mean(e[:, 0] + e[:, 1]) * 1e4)
    b3.plot(tt, em, color=c, label=rock)
b2.set_xlim(0, 360); b2.set_xticks(range(0, 361, 45))
b2.set_xlabel(r"$\theta$ [deg]"); b2.set_ylabel(r"$\varepsilon_{vol}$ [$10^{-2}$ %]")
b2.set_title("Volumetric strain on probe circle r = 0.065 m"); b2.legend(fontsize=7, ncol=2)
b3.axvspan(4000, 4050, color="0.85", zorder=0)
b3.set_xlabel("t [s]"); b3.set_ylabel(r"mean $\varepsilon_{vol}$ on r = 0.065 m [$10^{-2}$ %]")
b3.set_title("Mean volumetric strain on the probe circle"); b3.legend(fontsize=8)

for x in ax.flat:
    x.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(a.png, dpi=150)
print("saved", a.png)
