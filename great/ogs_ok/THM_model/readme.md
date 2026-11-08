# THM model – GREAT cell, Load A (07.11.2026)

This model combines the M model (`../M_model`, benchmark M1, Mollaali et al. 2026, JRMGE,
doi:10.1016/j.jrmge.2026.03.060) with the TH model (`../T_model/mesh/TH_granite_injection.prj`).
It runs in OGS 6.5.9 with THERMO_HYDRO_MECHANICS, 2D plane strain, monolithic Newton.

## Files
| File | Content |
|---|---|
| `make_great_cell_mesh.py` | GREAT cell mesh, **quadratic** triangles (T, p order 1; u order 2): rock r ≤ 0.097 m (MatID 0), sheath to 0.1 m (MatID 1), PEE/DSS arcs, probe circle r = 0.065 m, centre node; h = 1.5 mm at the centre, 6 mm at r = 0.065 m, 2 mm in the sheath (11889 nodes, 5804 elements) |
| `make_THM_prj.py` | writes `THM_greywacke_loadA.prj`, `THM_gneiss_loadA.prj` (`--load`, `--alpha_T`) |
| `plot_THM_overview.py` → `THM_loadA_overview.png` | centre T and p, radial T/p profiles, ε_vol on the probe circle, mean ε_vol vs. time |
| `THM_out_loadA.zip` | OGS results: domain at 5, 325, 2000, 4000, 4001, 4005, 4010, 4025, 4050, 4060, 4100, 4200, 4500, 5000 s; centre node and probe circle at every time step; logs |

## Schedule
- t = 0: mechanical Load A applied in the first step. Values from Table 4: PEE 1…8 = 10.0, 6.64, 4.46, 1.17, 1.0, 3.82, 7.80, 9.95 MPa; each DSS = mean of its neighbours
- 0–4000 s: outer boundary T rises linearly 293 → 363 K (0–325 s), then held at 363 K
- 4000–4050 s: injection at the centre node, 1 l/min = 1.667e-5 m³/s per m (as a mass rate of 1.667e-2 kg/s per m, because `equation_balance_type = mass`), T_inj = 293 K
- 4050–5000 s: shut-in
- Time steps: 5 s (0–400 s), 40 s (–4000 s), 1 s (injection), 2 s (–4150 s), 10 s (–5000 s); run time ≈ 13 min per rock

## Parameters
- Elasticity, Biot coefficient, porosity, solid density from Table 1:
  - Greywacke: 26.85 GPa / 0.27 / α = 0.8 / φ = 0.005
  - Gneiss: 83.9 GPa / 0.21 / α = 0.6 / φ = 0.001
  - Sheath: 0.1 GPa / 0.4 / α = 0
- Fluid: ρ = 1000 kg/m³ at p0, compressibility 4.4e-10 1/Pa (Table 3), μ = 1e-3 Pa s, c = 4180 J/(kg K), λ = 0.6 W/(m K)
- From the TH model:
  - k = 1.409e-14 m², p0 = 1 MPa, T0 = 293 K, drained outer boundary
  - Bulk λ = 2.8 W/(m K) and bulk ρc = 2700·790 J/(m³ K); the solid λ_s and c_s are back-calculated for each porosity
- Storage is not set directly. It follows from φ, α, the fluid compressibility and the elasticity.
- Initial effective stress α·p0, so that the initial state (p = p0, no load) has zero total stress

## Assumptions
- Linear thermal expansion coefficient of the solid 8e-6 1/K (`--alpha_T`); the fluid has no thermal expansion, as in the TH model
- Sheath hydraulically and thermally like the rock. This keeps the drained boundary of the TH model; the paper gives the sheath k = 1e-19 m².
- PEE arc 18°, DSS arc 4.5° (estimated from Fig. 4b)

## Checks
- Centre temperature vs. TH model: max. difference 0.42 K, caused by the coarser time steps during heating (40 s instead of 5 s)
- Injection: dp(r = 0.065 m) = 0.0812 MPa vs. Thiem 0.0811 MPa; dp at the centre node 1.12 MPa (TH model 1.09 MPa on a coarser centre mesh)
- ε_vol on the probe circle at t = 5 s vs. M1 model: max. difference < 0.01·10⁻² %
- Thermal strain at 4000 s: mean ΔT on the probe circle = 67.0 K, so 2(1+ν)α_T·ΔT + ε_mech = 11.15 (Greywacke) / 12.02 (Gneiss); THM gives 11.15 / 12.03 (units 10⁻² %)
