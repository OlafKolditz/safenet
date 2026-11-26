# Model description – UFZ, Scenario II, THM2 (Load A)

GREAT Cell Modelling Specification v2 (22.11.2026), Section 8 ("Model description") and Section 8.2.

| Item | UFZ model |
|---|---|
| Team | UFZ (Helmholtz Centre for Environmental Research, Leipzig) |
| Code, version | OpenGeoSys (OGS) 6.5.9 |
| Process formulation | `THERMO_HYDRO_MECHANICS`: quasi-static, fully saturated, single-phase thermo-poroelasticity, Biot effective stress, local thermal equilibrium, `equation_balance_type = mass`, no gravity |
| Coupling, solver | Monolithic Newton–Raphson; `PerComponentDeltaX` (NORM2) with absolute tolerances 1e-3 K, 10 Pa, 1e-8 m, 1e-8 m; Eigen SparseLU with scaling |
| Dimension | 2D plane strain, unit thickness 1 m; ε_vol = ε_xx + ε_yy |
| Elements | Quadratic triangles (Taylor–Hood): u order 2; p, T order 1; integration order 4 |
| Mesh | 11 889 nodes, 5 804 elements (gmsh); h = 1.5 mm at the centre, 6 mm at r_p = 0.065 m, 2 mm at 0.097 < r ≤ 0.1 m |
| Geometry | Disc R = 0.1 m, **no rubber sheath** (the mesh ring 0.097–0.1 m carries the rock properties); 16 PEE arcs (18°) and 16 DSS arcs (4.5°); θ ccw from +x, PEE 1 at 0° |
| Time stepping | Backward Euler, fixed steps: 10 s (0–2000 s), 1 s (2000–2120 s), 5 s (2120–5600 s); 1016 steps |
| Mechanical load | Load A, PEE 1…8 = 10.0, 6.64, 4.46, 1.17, 1.0, 3.82, 7.80, 9.95 MPa, DSS = mean of neighbours; linear ramp 0–1800 s, then held; rigid-body points u_x = u_y = 0 at (0, −0.09894), u_y = 0 at (−0.09894, 0) |
| Rock | Spec Tables 4/5: GW E = 26.85 GPa, ν = 0.27, α_B = 0.8, φ = 0.005, ρ_s = 2650, k = 2.58e-19 m², λ_s = 2.5, c_s = 850, α_s = 8e-6; GN E = 83.9 GPa, ν = 0.21, α_B = 0.6, φ = 0.001, ρ_s = 2750, k = 1e-19 m², λ_s = 3.0, c_s = 800, α_s = 7e-6; bulk λ by porosity mixing |
| Water | ρ(T), c(T): 4th-order polynomials fitted to IAPWS-IF97 (15–300 °C, p_abs ≈ 8.6 MPa; errors < 0.1 % / 0.4 %); μ(T) Vogel equation (error < 7 % at 300 °C); β_f = −(1/ρ) dρ/dT; λ_f = 0.6 W/(m K); compressibility 0; pressure dependence neglected |
| Initial conditions | T0 = 300 °C uniform, T_ref = T0; p = 0 (gauge); u = 0, no initial stress |
| Boundary conditions | Adiabatic perimeter; p = 0 on the perimeter; normal tractions on PEE/DSS arcs |
| Injection | Nodal mass source at the centre node: Q2D = Q0/B = 8.333e-9 m³/s per m (Q0 = 0.1 ml/min, B = 0.2 m), mass rate ρ_w Q2D with ρ_w = 1000 kg/m³; linear ramp 2000–2120 s, constant to 5600 s; no pressure limit reached/applied |
| **Inlet condition** | **Fixed temperature 15 °C (Dirichlet) at the centre node for t > 2000.5 s** (point source; not an advective enthalpy flux). This acts as a heat sink: cooling is conduction-controlled (cf. spec Sect. 7.2 and 10.2) |
| Run time | ≈ 52 min per case on one core |

## Output files (spec Table 16)
| File | Content |
|---|---|
| `UFZ_<case>_points.csv` | t, T_P0, p_P0, eps_vol_P0, …, T_P270, p_P270, eps_vol_P270; every time step 0–5600 s (1017 rows) |
| `UFZ_<case>_inlet.csv` | t, p_in, T_in, Q_ml_min, V_cum_ml at the centre node; every time step |
| `UFZ_<case>_profiles.csv` | t, theta_deg, T, p, eps_vol on r_p = 0.065 m, θ = 0…359° (1°), t = 2000, 2060, 3800, 5600 s |
| `UFZ_<case>_field_t<time>.vtu` | full field at 2000, 2060, 3800, 5600 s: u_magnitude [m], sigma_trace [MPa] (total stress σ = σ′ − α_B p I, 3D trace incl. σ_zz), sigma_eff_trace [MPa], eps_vol [–], p [MPa], q [m/s], T [°C], u [m] |
| `UFZ_THM2_scII_histories.png`, `UFZ_THM2_scII_profiles.png` | overview plots (`workflow/plot_THM2_scII.py`) |

Units: t in s, T in °C, p in MPa (gauge, p_abs = p + 8.6 MPa), strain as fraction (tension positive), Q in ml/min, V in ml.
Probe-circle values: nodal values at the element corner nodes (every 5°), linear interpolation in θ for the 1° profiles;
the same rule for T, p and ε_vol. P0, P90, P180, P270 are nodes. `events` files are not applicable (no fracture).

## Key results (t = 5600 s)
| Case | p_in [MPa] | min p_in [MPa] | T_P0 [°C] | p_P0 [MPa] | ε_vol,P0 [µε] | ε_vol,P90 [µε] |
|---|---|---|---|---|---|---|
| THM2-GW-A-T300 | 7.46 | −3.01 | 263.5 | 0.215 | −1086 | −930 |
| THM2-GN-A-T300 | 19.93 | −4.48 | 254.8 | 0.682 | −912 | −852 |

- Cumulative injected volume 5.90 ml (nominal 6.0 ml minus 0.1 ml for the 120 s ramp).
- Symmetry: P0 = P180 and P90 = P270 within 0.007 K and 0.1 µε.
- Effect of removing the sheath (vs. OGS model v2 with sheath): T and p unchanged (< 0.1 %); ε_vol at P0 9 µε (GW) and 6 µε (GN) less compressive at 5600 s (< 1 %).
- p_in is the point-source value at the centre node and depends on the element size there (1.5 mm).
- The short negative inlet pressure pulse at 2001–2032 s (GW) and 2001–2018 s (GN) (cooling contraction of the pore water at the inlet) brings p below 0 gauge, i.e. below the saturation pressure at 300 °C; the single-phase model ignores boiling (spec Sect. 10.2).

## Deviations from the specification
- Isothermal HM controls (HM2-GW-A-T15, HM2-GN-A-T15) are not part of this delivery.
- PEE/DSS arc widths 18°/4.5° are estimates (spec open item).
- Temperature-dependent water properties use fits, not full IAPWS-IF97; pressure dependence neglected.
