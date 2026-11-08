"""
Figures for the GREAT cell THM model report: geometry, mesh, boundary conditions.
Reads the OGS meshes of THM_model (domain.vtu, PEE_i.vtu, DSS_i.vtu, ...).

Usage: python make_report_figures.py [--mesh THM_model] [--out figures]
"""
import argparse
import os

import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Rectangle

R_OUT, R_ROCK, R_PROBE, R_FIX = 0.100, 0.097, 0.065, 0.09894
LOAD_A = [10.0, 6.64, 4.46, 1.17, 1.0, 3.82, 7.80, 9.95]       # MPa, PEE 1..8 (& 1a..8a)
DSS_A = [(LOAD_A[i] + LOAD_A[(i + 1) % 8]) / 2 for i in range(8)]
C_ROCK, C_SHEATH = "#d9d4c7", "#7a9cc6"
C_PEE, C_DSS = "#c0392b", "#2c7fb8"

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "savefig.dpi": 300,
                     "font.family": "serif", "mathtext.fontset": "cm"})


def load(mdir, name):
    return meshio.read(os.path.join(mdir, name))


def segs(m):
    """line (or line3) cells -> list of polylines in (x, y)"""
    c = m.cells[0].data
    order = [0, 2, 1] if c.shape[1] == 3 else [0, 1]
    return [m.points[row[order], :2] for row in c]


def tri_corners(m):
    return m.points[m.cells[0].data[:, :3], :2]


def arc_mid(m):
    p = m.points[:, :2]
    ang = np.arctan2(p[:, 1], p[:, 0])
    # split into the two opposite arcs (i and ia) by sign of the projection on the mean direction
    z = np.exp(1j * ang)
    k = np.argmax(np.abs(z))
    d = z[k]
    a1 = ang[np.real(z * np.conj(d)) > 0]
    a2 = ang[np.real(z * np.conj(d)) <= 0]
    mid = lambda a: np.angle(np.mean(np.exp(1j * a)))
    return [mid(a1), mid(a2)]


def circle(ax, r, **kw):
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(r * np.cos(t), r * np.sin(t), **kw)


def fig_geometry(mdir, out):
    dom = load(mdir, "domain.vtu")
    mat = dom.cell_data["MaterialIDs"][0]
    tri = tri_corners(dom)
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={"width_ratios": [1.15, 1]})
    for k, ax in enumerate(axs):
        ax.add_collection(PolyCollection(tri[mat == 0], facecolors=C_ROCK, edgecolors="none"))
        ax.add_collection(PolyCollection(tri[mat == 1], facecolors=C_SHEATH, edgecolors="none"))
        for i in range(1, 9):
            for kind, col in (("PEE", C_PEE), ("DSS", C_DSS)):
                m = load(mdir, f"{kind}_{i}.vtu")
                ax.add_collection(LineCollection(segs(m), colors=col, linewidths=3.0 if k else 2.2,
                                                 capstyle="butt"))
        circle(ax, R_PROBE, color="k", lw=0.8, ls="--")
        ax.plot(0, 0, "k+", ms=8, mew=1.2)
        for fx, fy, lab in ((0, -R_FIX, r"$u_x=u_y=0$"), (-R_FIX, 0, r"$u_y=0$")):
            ax.plot(fx, fy, marker="^" if fy else ">", color="k", ms=5 if k == 0 else 8)
        ax.set_aspect("equal")
        ax.set_xlabel("$x$ (m)")
    ax = axs[0]
    for i in range(1, 9):
        m1, m2 = arc_mid(load(mdir, f"PEE_{i}.vtu"))
        for a, s in ((m1, ""), (m2, "a")):
            ax.text(1.13 * R_OUT * np.cos(a), 1.13 * R_OUT * np.sin(a), f"{i}{s}",
                    ha="center", va="center", fontsize=7.5, color=C_PEE)
    a_p = np.radians(140)
    ax.annotate("probe circle $r=0.065$ m", (R_PROBE * np.cos(a_p), R_PROBE * np.sin(a_p)),
                (-0.005, 0.022), fontsize=7.5, ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.annotate("centre node\n(injection)", (0, 0), (0.012, -0.03), fontsize=7.5,
                arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.text(0.0, 0.074, "rock $r\\leq 0.097$ m", ha="center", fontsize=7.5)
    lim = 0.122
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_ylabel("$y$ (m)")
    ax.add_patch(Rectangle((-0.1035, -0.024), 0.024, 0.048, fill=False, lw=0.6, ls=":"))
    ax.set_title("(a) model domain, PEE numbering", fontsize=9)
    ax = axs[1]
    ax.set_xlim(-0.1035, -0.0795); ax.set_ylim(-0.024, 0.024)
    ax.text(-0.0835, -0.012, "rock (MatID 0)", rotation=90, ha="center", va="center", fontsize=8)
    ax.annotate("sheath\n(MatID 1)", (-0.0985, -0.008), (-0.0945, -0.008), fontsize=7.5,
                va="center", arrowprops=dict(arrowstyle="-", lw=0.5))
    a_d = np.radians(180 - 9 - 2.25)
    ax.annotate("DSS 8 (4.5°)", (R_OUT * np.cos(a_d), R_OUT * np.sin(a_d)), (-0.094, 0.0205),
                fontsize=7.5, color=C_DSS, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=C_DSS))
    ax.annotate("PEE 1a (18°)", (-0.1, 0.008), (-0.094, 0.011), fontsize=7.5, color=C_PEE,
                va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color=C_PEE))
    ax.annotate(r"$P_\mathrm{left}$: $u_y=0$", (-R_FIX, 0), (-0.094, 0.003), fontsize=7.5,
                va="center", arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.set_title("(b) detail at PEE 1a (sheath 3 mm)", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_geometry.pdf")); fig.savefig(os.path.join(out, "fig_geometry.png"))
    plt.close(fig)


def fig_mesh(mdir, out):
    dom = load(mdir, "domain.vtu")
    mat = dom.cell_data["MaterialIDs"][0]
    tri = tri_corners(dom)
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw={"width_ratios": [1.3, 1, 1]})
    views = [(-0.105, 0.105, -0.105, 0.105, "(a) full mesh"),
             (-0.012, 0.012, -0.012, 0.012, "(b) centre, $h=1.5$ mm"),
             (0.080, 0.101, -0.012, 0.012, "(c) rock–sheath, $h=2$ mm")]
    for ax, (x0, x1, y0, y1, title) in zip(axs, views):
        lw = 0.1 if title.startswith("(a)") else 0.35
        ax.add_collection(PolyCollection(tri[mat == 0], facecolors=C_ROCK, edgecolors="k",
                                         linewidths=lw))
        ax.add_collection(PolyCollection(tri[mat == 1], facecolors=C_SHEATH, edgecolors="k",
                                         linewidths=lw))
        if not title.startswith("(a)"):
            circle(ax, R_PROBE, color="r", lw=0.8)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
        ax.set_title(title, fontsize=9)
        ax.tick_params(labelsize=7)
        if not title.startswith("(a)"):
            ax.set_xticks(np.round(np.linspace(x0, x1, 3)[:], 3) if x0 < 0 else [0.08, 0.09, 0.10])
        ax.set_xlabel("$x$ (m)")
    axs[0].set_ylabel("$y$ (m)")
    for (x0, x1, y0, y1, _) in views[1:]:
        axs[0].add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="r", lw=0.8))
    axs[1].plot(0, 0, "r+", ms=8)
    n_nodes = len(dom.points); n_el = len(tri)
    fig.suptitle(f"{n_nodes} nodes, {n_el} quadratic triangles (corner edges shown)",
                 fontsize=8, y=0.02, va="bottom")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(os.path.join(out, "fig_mesh.pdf")); fig.savefig(os.path.join(out, "fig_mesh.png"))
    plt.close(fig)
    return n_nodes, n_el


def fig_bc(mdir, out):
    fig = plt.figure(figsize=(7.2, 3.5))
    ax = fig.add_axes([0.02, 0.08, 0.44, 0.84])
    ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(plt.Circle((0, 0), R_OUT, fc=C_ROCK, ec="k", lw=0.6))
    circle(ax, R_PROBE, color="k", lw=0.6, ls="--")
    smax = max(LOAD_A)
    for kind, vals, col in (("PEE", LOAD_A, C_PEE), ("DSS", DSS_A, C_DSS)):
        for i, v in enumerate(vals, 1):
            for a in arc_mid(load(mdir, f"{kind}_{i}.vtu")):
                L = 0.012 + 0.045 * v / smax
                x1, y1 = R_OUT * np.cos(a), R_OUT * np.sin(a)
                ax.annotate("", (x1, y1), (x1 + L * np.cos(a), y1 + L * np.sin(a)),
                            arrowprops=dict(arrowstyle="-|>", color=col,
                                            lw=1.3 if kind == "PEE" else 0.7,
                                            mutation_scale=7 if kind == "PEE" else 5))
                if kind == "PEE":
                    rr = R_OUT + L + 0.012
                    ax.text(rr * np.cos(a), rr * np.sin(a), f"{v:g}", fontsize=7,
                            ha="center", va="center", color=col)
    ax.plot(0, 0, "o", color="#1a9850", ms=5)
    ax.text(0.004, -0.014, "$Q_\\mathrm{inj}$, $T_\\mathrm{inj}$", fontsize=7.5, color="#1a9850")
    lim = 0.168
    ax.text(-lim, -lim + 0.012, "outer boundary: $T=T_\\mathrm{outer}(t)$, $p=p_0=1$ MPa",
            fontsize=7, va="bottom")
    lim = 0.168
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_title("(a) Load A: normal traction (MPa)", fontsize=9)
    ax.text(-lim, -lim, "PEE (red), DSS = mean of neighbours (blue)", fontsize=7, va="bottom")

    t = np.array([0, 325, 4000, 5000])
    T = np.array([293, 363, 363, 363])
    ax1 = fig.add_axes([0.56, 0.17, 0.34, 0.72])
    ax1.plot(t, T, color=C_PEE, lw=1.4)
    ax1.set_xlabel("time $t$ (s)"); ax1.set_ylabel("$T_\\mathrm{outer}$ (K)", color=C_PEE)
    ax1.set_xlim(0, 5000); ax1.set_ylim(285, 375)
    ax2 = ax1.twinx()
    tq = [0, 4000, 4000, 4050, 4050, 5000]
    q = [0, 0, 1, 1, 0, 0]
    ax2.fill_between(tq, q, step=None, color="#1a9850", alpha=0.35, lw=0)
    ax2.plot(tq, q, color="#1a9850", lw=1.0)
    ax2.set_ylabel("$Q_\\mathrm{inj}$ (l/min)", color="#1a9850"); ax2.set_ylim(0, 1.25)
    for x0, x1, lab in ((0, 4000, "heating"), (4000, 4050, ""), (4050, 5000, "shut-in")):
        ax1.axvspan(x0, x1, color="0.5", alpha=0.06 if lab != "" else 0.0)
        if lab:
            ax1.text(0.5 * (x0 + x1), 368, lab, ha="center", fontsize=7.5)
    ax1.annotate("injection\n4000–4050 s", (4025, 300), (2600, 300), fontsize=7.5,
                 arrowprops=dict(arrowstyle="-|>", lw=0.6))
    ax1.axvline(0, color="k", lw=2)
    ax1.text(80, 298, "Load A\nat $t=0$", fontsize=7.5)
    ax1.set_title("(b) thermal and hydraulic schedule", fontsize=9)
    fig.savefig(os.path.join(out, "fig_bc.pdf")); fig.savefig(os.path.join(out, "fig_bc.png"))
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", default="THM_model")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fig_geometry(a.mesh, a.out)
    print(fig_mesh(a.mesh, a.out))
    fig_bc(a.mesh, a.out)
