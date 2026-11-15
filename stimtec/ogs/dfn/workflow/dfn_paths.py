"""dfn_paths.py  -  central folder layout of the STIMTEC DFN model (dfn2)

    dfn2/
        simulation/   OGS project (*.prj) and meshes (*.vtu)
        results/      OGS output (*.pvd, *.vtu, logs)
        workflow/     Python tools (this folder), calibration logs, plots

Same layout as the single-fracture models in BH10_20180716_40.6.
All paths are built relative to this file, so the tools work from any working directory.
"""
import os

WORKFLOW_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.path.normpath(os.path.join(WORKFLOW_DIR, ".."))       # dfn2
STIMTEC_DIR = os.path.normpath(os.path.join(MODEL_DIR, ".."))        # stimtec


def default_results_dir(script_dir=WORKFLOW_DIR):
    """Results folder: ../results relative to the workflow folder."""
    return os.path.normpath(os.path.join(script_dir, "..", "results"))


def default_simulation_dir(script_dir=WORKFLOW_DIR):
    """OGS project folder: ../simulation relative to the workflow folder."""
    return os.path.normpath(os.path.join(script_dir, "..", "simulation"))


def bh10_40_file(name):
    """Data file of pumping test BH10_20180716_40.6 (measured pressure, injection stages).
    Looked up in stimtec/BH10_20180716_40.6/workflow, fallback: this workflow folder."""
    for d in (os.path.join(STIMTEC_DIR, "BH10_20180716_40.6", "workflow"), WORKFLOW_DIR):
        f = os.path.join(d, name)
        if os.path.isfile(f):
            return f
    raise FileNotFoundError(f"{name} neither in BH10_20180716_40.6/workflow nor in {WORKFLOW_DIR}")


RESULTS_DIR = default_results_dir()
SIMULATION_DIR = default_simulation_dir()
MEAS_NAME = "bh10_40_gemessen_digitalisiert.csv"
STEPS_NAME = "injektion_stufen_40.csv"
