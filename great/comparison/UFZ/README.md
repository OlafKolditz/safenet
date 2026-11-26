# UFZ contribution – GREAT cell Scenario II, THM2 (Load A)

Team: UFZ (Helmholtz Centre for Environmental Research), contact O. Kolditz
Specification: GREAT Cell Modelling Specification v2 (22.11.2026), Scenario II, outputs per Section 8.2
Cases: `THM2-GW-A-T300` (Greywacke, Phase 1 mandatory) and `THM2-GN-A-T300` (Gneiss, optional)

## Folders
| Folder | Content |
|---|---|
| `simulation/` | OGS project files (`THM2-GW-A-T300.prj`, `THM2-GN-A-T300.prj`), meshes (`domain.vtu` + 21 boundary/point meshes, `great_cell.msh`), generators `make_THM2_scII_prj.py` and `make_great_cell_mesh.py`, OGS logs `log_<case>.txt` |
| `results/` | Output files of Section 8.2 (Table 16) for both cases, model description `model_description.md`, overview plots |
| `workflow/` | `export_THM2_scII.py` (OGS output → spec files), `plot_THM2_scII.py` (plots from the CSV files) |

## Workflow
```
cd simulation
python make_great_cell_mesh.py                  # optional: meshes are included
python make_THM2_scII_prj.py --load A           # writes THM2-GW-A-T300.prj, THM2-GN-A-T300.prj
ogs -o out THM2-GW-A-T300.prj                   # OGS 6.5.9, about 52 min each on one core
ogs -o out THM2-GN-A-T300.prj
cd ../workflow
python export_THM2_scII.py ../simulation/out ../results
python plot_THM2_scII.py ../results --ref ../../../THM_model_v2/results
```
Python packages: numpy, pandas, matplotlib, pyvista (gmsh for the mesh).

## Not included
- Isothermal HM controls (`HM2-GW-A-T15`, `HM2-GN-A-T15`) required by the specification: not rerun in this set.
  They can be generated with `python make_THM2_scII_prj.py --T0 15`.
- Raw OGS output (`simulation/out/`, about 76 MB for both cases), reproducible from the project files.
