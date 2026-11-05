import colorsys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot

_PLOT_RC = {
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "legend.fontsize": 8,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "lines.linewidth": 1.5,
    "axes.linewidth": 0.8,
    "grid.linewidth": 0.4,
    "savefig.dpi": 600,
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

_h, _lum, _s = colorsys.rgb_to_hls(*[int(c, 16) / 255 for c in ("c0", "39", "2b")])
_lum = max(0.2, min(0.9, 0.7 * (_lum or 0.5)))
_COLOR = "#{:02x}{:02x}{:02x}".format(*[int(x * 255) for x in colorsys.hls_to_rgb(_h, _lum, _s)])


def plot_vol_strain(mesh, label, out_dir, probe_radius=0.065, ylim=(-7.5, 2.5)):
    out_dir = Path(out_dir)
    radii = np.linalg.norm(mesh.cell_centers().points[:, :2], axis=1)
    inner = mesh.extract_cells(radii < probe_radius)
    edge  = inner.extract_feature_edges()

    eps_v   = ot.variables.strain.trace.transform(edge)
    phi     = np.arctan2(edge.points[:, 1], edge.points[:, 0])
    phi[phi < 0] += 2 * np.pi
    idx     = np.argsort(phi)
    phi_deg = np.degrees(phi[idx])
    eps_pct = np.array(eps_v)[idx] * 100

    plt.rcParams.update(_PLOT_RC)
    mk = max(1, len(phi_deg) // 10)

    fig, ax = plt.subplots(1, 1, figsize=(7.0, 4.0))
    ax.set_title("Load case: A", fontsize=10, pad=10)
    ax.plot(phi_deg, eps_pct, label=label,
            color=_COLOR, linestyle="-", marker="^", markersize=6,
            markevery=mk, markerfacecolor="none", markeredgewidth=1.2)
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.set_ylim(*ylim)
    ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.7)
    ax.set_ylabel(r"Volumetric strain $\varepsilon_{\mathrm{vol}}$ / %", labelpad=5)
    ax.set_xlabel(r"Angle $\theta$ / $^\circ$", labelpad=5)
    ax.minorticks_on()
    ax.tick_params(which="minor", length=3, width=0.6)
    ax.tick_params(which="major", length=5, width=0.8)
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax.legend(loc="lower right", bbox_to_anchor=(0.98, 0.098),
              frameon=True, framealpha=0.95, edgecolor="0.3", borderpad=0.6)

    fig.tight_layout()
    tag   = label.split(",")[1].strip().replace(" ", "_") if "," in label else label.replace(" ", "_")
    fname = out_dir / f"volume_strain_plot_{tag}_subplots.pdf"
    fig.savefig(fname, dpi=600, bbox_inches="tight", pad_inches=0.05)
    plt.show()


if __name__ == "__main__":
    import sys
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../_out")

    cases = [
        ("M1_SD_A_Greywacke.pvd",    "Greywacke, SD, Load A"),
        ("HM1_HM_A_Greywacke.pvd",   "Greywacke, HM, Load A"),
        ("THM1_THM_A_Greywacke.pvd", "Greywacke, THM, Load A"),
    ]
    for pvd, label in cases:
        p = out_dir / pvd
        if p.exists():
            ms = ot.MeshSeries(str(p))
            plot_vol_strain(ms[-1], label, out_dir)
        else:
            print(f"skipping {pvd} (not found in {out_dir})")
