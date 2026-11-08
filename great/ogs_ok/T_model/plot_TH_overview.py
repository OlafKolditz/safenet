"""Overview plot for TH_granite_injection.prj (OGS 6.5.9).

Usage: python plot_TH_overview.py <output_dir> <png>
Panels: (a) centre node T and dp over full run, (b) zoom 3950-5000 s,
        (c) radial temperature profiles, (d) radial pressure profiles (+ Thiem).
"""
import sys, glob, re
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = sys.argv[1] if len(sys.argv) > 1 else "out"
png = sys.argv[2] if len(sys.argv) > 2 else "TH_granite_injection_overview.png"

p0 = 1e6
R, r_in = 0.1, 0.066
Q, mu, k = 1.6666667e-5, 1e-3, 1.409e-14


def files(mesh):
    fs = glob.glob(f"{out}/TH_granite_{mesh}_ts_*_t_*.vtu")
    t = [float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)) for f in fs]
    o = np.argsort(t)
    return np.array(t)[o], [fs[i] for i in o]


# centre node time series
tc, fc = files("center")
Tc, pc = [], []
for f in fc:
    m = pv.read(f)
    Tc.append(m.point_data["temperature"][0])
    pc.append(m.point_data["pressure"][0])
Tc, dpc = np.array(Tc), (np.array(pc) - p0) / 1e6

# domain profiles
td, fd = files("domain")
dom = {round(t, 3): f for t, f in zip(td, fd)}


def profile(t):
    m = pv.read(dom[t])
    r = np.hypot(m.points[:, 0], m.points[:, 1])
    o = np.argsort(r)
    return r[o], m.point_data["temperature"][o], (m.point_data["pressure"][o] - p0) / 1e6


T_times = [4000, 4010, 4025, 4050, 4100, 4200, 4500, 5000]
p_times = [4000, 4001, 4010, 4025, 4050, 4060, 4100]
cmap = plt.get_cmap("tab10")
col = {t: cmap(i) for i, t in enumerate(T_times)}          # same colours as before
pcol = {t: plt.get_cmap("viridis")(i / (len(p_times) - 1)) for i, t in enumerate(p_times)}

fig, ax = plt.subplots(2, 2, figsize=(14, 9.5))
(a, b), (c, d) = ax


def ts_panel(axx, xlim, title):
    axx.plot(tc, Tc, color="tab:red", label="T")
    axx.set_ylabel("T [K]")
    axx.set_xlabel("t [s]")
    axx.axvspan(4000, 4050, color="0.85", zorder=0)
    ax2 = axx.twinx()
    ax2.plot(tc, dpc, color="tab:blue")
    ax2.set_ylabel("p - p0 [MPa]", color="tab:blue")
    axx.set_xlim(*xlim)
    axx.set_title(title)


ts_panel(a, (-200, 5200), "Centre node: temperature (red) and pressure (blue)")
ts_panel(b, (3950, 5050), "Centre node, zoom 3950-5000 s")

for t in T_times:
    r, T, _ = profile(t)
    c.plot(r, T, ".", ms=1.5, color=col[t], label=f"t = {t} s")
c.axvline(r_in, ls="--", lw=0.8, color="0.4")
c.set_xlabel("r [m]"); c.set_ylabel("T [K]")
c.set_title("Radial temperature profiles")
c.legend(markerscale=6, fontsize=8, loc="lower right")

# draw large-to-small markers so overlapping (quasi-steady) profiles stay visible
for i, t in enumerate(p_times):
    r, _, dp = profile(t)
    d.plot(r[r > 0], dp[r > 0], "o", ms=7 - 0.8 * i, mfc="none" if i % 2 else pcol[t],
           mec=pcol[t], mew=0.8, ls="none", label=f"t = {t} s")
d.set_xscale("log")
d.text(0.03, 0.06, "4001-4050 s: identical profiles (pressure diffusivity\n"
       "k/(mu S) ~ 1 m2/s -> quasi-steady within < 0.1 s)\n"
       "4000 s and after shut-in (>= 4060 s): p = p0",
       transform=d.transAxes, fontsize=8, va="bottom")
rr = np.logspace(-3, np.log10(R), 400)
d.plot(rr, Q * mu / (2 * np.pi * k) * np.log(R / rr) / 1e6, "k-", lw=0.8,
       label="Thiem (steady)")
d.axvline(r_in, ls="--", lw=0.8, color="0.4")
d.set_xlabel("r [m]"); d.set_ylabel("p - p0 [MPa]")
d.set_title("Radial pressure profiles")
d.set_ylim(-0.05, 1.15)
d.legend(markerscale=1, fontsize=8, loc="upper right")

fig.tight_layout()
fig.savefig(png, dpi=150)
print("saved", png)
