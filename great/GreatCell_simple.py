#!/Volumes/OpenGeoSys/build/release/.venv/bin/python
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: OGS Build (.venv)
#     language: python
#     name: ogs-build-env
# ---

# %% [raw] magic_args="[raw]"
# +++
# title = "GREAT cell — M1, HM1, THM1 (intact sample), Load A"
# date = "2026-11-05"
# author = "Mostafa Mollaali, Olaf Kolditz"
# web_subsection = "small-deformations"
# +++

# %% [markdown]
# # GREAT cell — M1, HM1, THM1 — Load A
#
# Three simulations of the intact Greywacke sample under load case A:
#
# - **M1**   — small-deformation mechanics only
# - **HM1**  — hydro-mechanics with fluid injection at the centre (Inlet)
# - **THM1** — thermo-hydro-mechanics: cold water (20 °C) injected into hot rock (100 °C)
#
# All share the same mesh (`mesh/`).
# Post-processing functions live in `postprocess.py`.

# %% [markdown]
# ## Setup

# %%
import os
from pathlib import Path

import matplotlib as mpl
import ogstools as ot
from postprocessing.snapshots import plot_snapshots
from postprocessing.vol_strain import plot_vol_strain

mpl.rcdefaults()
mpl.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "axes.labelsize": 26,
        "axes.titlesize": 24,
        "legend.fontsize": 18,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
    }
)

# %%
OGS_BIN = None
out_dir  = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# Set any flag to False to skip that simulation and use existing _out/ results
run_M1   = True
run_HM1  = True
run_THM1 = True

# %% [markdown]
# ---
# # M1 — Small-deformation mechanics

# %% [markdown]
# ## Run OGS

# %%
if run_M1:
    prj = ot.Project(input_file="M1_SD_A.prj", output_dir=str(out_dir))
    prj.run_model(logfile=out_dir / "out_M1_SD_A.log", path=OGS_BIN, args="-m mesh")

# %% [markdown]
# ## Post-processing

# %%
ms_m1   = ot.MeshSeries(str(out_dir / "M1_SD_A_Greywacke.pvd"))
mesh_m1 = ms_m1[-1]
print(f"M1  — t = {ms_m1.timevalues[-1]:.0f} s  |  fields: {list(mesh_m1.point_data.keys())}")

# %% [markdown]
# ### Volumetric strain, circle probe

# %%
plot_vol_strain(mesh_m1, "Greywacke, SD, Load A", out_dir)

# %% [markdown]
# ### Field snapshots

# %%
displacement = ot.variables.displacement.replace(output_unit="mm")

plot_snapshots(
    mesh_m1,
    fields=[
        (displacement,               r"$\mathbf{u}$ / mm",           {"vmin": 0,     "vmax": 0.25}),
        (ot.variables.stress.trace,  r"$\mathrm{tr}(\sigma)$ / MPa", {"vmin": -20,   "vmax": 10}),
        (ot.variables.strain.trace,  r"$\mathrm{tr}(\varepsilon)$",  {"vmin": -0.05, "vmax": 0}),
    ],
    out_dir=out_dir,
    filename="M1_SD_A_field_snapshots.pdf",
)

# %% [markdown]
# ---
# # HM1 — Hydro-mechanics with fluid injection

# %% [markdown]
# ## Run OGS

# %%
if run_HM1:
    prj = ot.Project(input_file="HM1_HM_A.prj", output_dir=str(out_dir))
    prj.run_model(logfile=out_dir / "out_HM1_HM_A.log", path=OGS_BIN, args="-m mesh")

# %% [markdown]
# ## Post-processing

# %%
ms_hm1   = ot.MeshSeries(str(out_dir / "HM1_HM_A_Greywacke.pvd"))
mesh_hm1 = ms_hm1[-1]
print(f"HM1 — t = {ms_hm1.timevalues[-1]:.1f} s  |  fields: {list(mesh_hm1.point_data.keys())}")

# %% [markdown]
# ### Volumetric strain, circle probe

# %%
plot_vol_strain(mesh_hm1, "Greywacke, HM, Load A", out_dir)

# %% [markdown]
# ### Field snapshots

# %%
plot_snapshots(
    mesh_hm1,
    fields=[
        (displacement,               r"$\mathbf{u}$ / mm",           {"vmin": 0,     "vmax": 0.25}),
        (ot.variables.pressure,      r"$p$ / MPa",                   {"vmin": 0.1,   "vmax": 5}),
        (ot.variables.stress.trace,  r"$\mathrm{tr}(\sigma)$ / MPa", {"vmin": -20,   "vmax": 10}),
        (ot.variables.strain.trace,  r"$\mathrm{tr}(\varepsilon)$",  {"vmin": -0.05, "vmax": 0}),
    ],
    out_dir=out_dir,
    filename="HM1_HM_A_field_snapshots.pdf",
)

# %% [markdown]
# ---
# # THM1 — Thermo-hydro-mechanics: cold water injection into hot rock
#
# Rock initially at **100 °C** (373.15 K), cold water injected at **20 °C** (293.15 K).
# p_ic = 1 MPa keeps water liquid. Thermal properties from literature (Vosteen & Schellschmidt 2003).

# %% [markdown]
# ## Run OGS

# %%
if run_THM1:
    prj = ot.Project(input_file="THM1_THM_A.prj", output_dir=str(out_dir))
    prj.run_model(logfile=out_dir / "out_THM1_THM_A.log", path=OGS_BIN, args="-m mesh")

# %% [markdown]
# ## Post-processing

# %%
ms_thm1   = ot.MeshSeries(str(out_dir / "THM1_THM_A_Greywacke.pvd"))
mesh_thm1 = ms_thm1[-1]
print(f"THM1 — t = {ms_thm1.timevalues[-1]:.1f} s  |  fields: {list(mesh_thm1.point_data.keys())}")

# %% [markdown]
# ### Volumetric strain, circle probe

# %%
plot_vol_strain(mesh_thm1, "Greywacke, THM, Load A", out_dir)

# %% [markdown]
# ### Field snapshots

# %%
plot_snapshots(
    mesh_thm1,
    fields=[
        (displacement,                r"$\mathbf{u}$ / mm",           {"vmin": 0,     "vmax": 0.25}),
        (ot.variables.pressure,       r"$p$ / MPa",                   {"vmin": 1.0,   "vmax": 5}),
        (ot.variables.temperature,    r"$T$ / °C",                    {"vmin": 20,    "vmax": 100}),
        (ot.variables.stress.trace,   r"$\mathrm{tr}(\sigma)$ / MPa", {"vmin": -20,   "vmax": 10}),
        (ot.variables.strain.trace,   r"$\mathrm{tr}(\varepsilon)$",  {"vmin": -0.05, "vmax": 0}),
    ],
    out_dir=out_dir,
    filename="THM1_THM_A_field_snapshots.pdf",
)
