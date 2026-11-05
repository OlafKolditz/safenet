from pathlib import Path

import matplotlib.pyplot as plt
import ogstools as ot

_DOMAIN_LIMIT = 0.098


def plot_snapshots(mesh, fields, out_dir, filename):
    out_dir = Path(out_dir)
    n = len(fields)
    fig, axs = plt.subplots(1, n, figsize=(12 * n, 10))
    if n == 1:
        axs = [axs]

    for ax, (var, title, opts) in zip(axs, fields, strict=True):
        ot.plot.contourf(mesh, var, fig=fig, ax=ax, **opts)
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xlabel("x / m", labelpad=4)
        ax.set_ylabel("y / m", labelpad=4)
        ax.set_xlim(-_DOMAIN_LIMIT, _DOMAIN_LIMIT)
        ax.set_ylim(-_DOMAIN_LIMIT, _DOMAIN_LIMIT)

    fig.tight_layout()
    fig.savefig(out_dir / filename, dpi=600, bbox_inches="tight", pad_inches=0.05)
    plt.show()


if __name__ == "__main__":
    import sys
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../_out")

    displacement = ot.variables.displacement.replace(output_unit="mm")

    cases = [
        (
            "M1_SD_A_Greywacke.pvd",
            "M1_SD_A_field_snapshots.pdf",
            [
                (displacement,               r"$\mathbf{u}$ / mm",           {"vmin": 0,     "vmax": 0.25}),
                (ot.variables.stress.trace,  r"$\mathrm{tr}(\sigma)$ / MPa", {"vmin": -20,   "vmax": 10}),
                (ot.variables.strain.trace,  r"$\mathrm{tr}(\varepsilon)$",  {"vmin": -0.05, "vmax": 0}),
            ],
        ),
        (
            "HM1_HM_A_Greywacke.pvd",
            "HM1_HM_A_field_snapshots.pdf",
            [
                (displacement,               r"$\mathbf{u}$ / mm",           {"vmin": 0,     "vmax": 0.25}),
                (ot.variables.pressure,      r"$p$ / MPa",                   {"vmin": 0.1,   "vmax": 5}),
                (ot.variables.stress.trace,  r"$\mathrm{tr}(\sigma)$ / MPa", {"vmin": -20,   "vmax": 10}),
                (ot.variables.strain.trace,  r"$\mathrm{tr}(\varepsilon)$",  {"vmin": -0.05, "vmax": 0}),
            ],
        ),
        (
            "THM1_THM_A_Greywacke.pvd",
            "THM1_THM_A_field_snapshots.pdf",
            [
                (displacement,               r"$\mathbf{u}$ / mm",           {"vmin": 0,     "vmax": 0.25}),
                (ot.variables.pressure,      r"$p$ / MPa",                   {"vmin": 1.0,   "vmax": 8}),
                (ot.variables.temperature,   r"$T$ / °C",                    {"vmin": 20,    "vmax": 100}),
                (ot.variables.stress.trace,  r"$\mathrm{tr}(\sigma)$ / MPa", {"vmin": -20,   "vmax": 10}),
                (ot.variables.strain.trace,  r"$\mathrm{tr}(\varepsilon)$",  {"vmin": -0.05, "vmax": 0}),
            ],
        ),
    ]

    for pvd, filename, fields in cases:
        p = out_dir / pvd
        if p.exists():
            ms = ot.MeshSeries(str(p))
            plot_snapshots(ms[-1], fields, out_dir, filename)
        else:
            print(f"skipping {pvd} (not found in {out_dir})")
