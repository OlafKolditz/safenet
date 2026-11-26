"""
OGS THM model (THERMO_HYDRO_MECHANICS, 2D plane strain) of the GREAT cell sample - Scenario II.
Comparison set for the GREAT Cell Modelling Specification v2 (22.11.2026), Scenario II, Section 8.2.
Derived from ../../THM_model_v2/make_THM_prj_v2.py; differences to model v2:
  - no rubber sheath (spec v2, Sect. 1.2/7.2): the mesh ring 0.097 < r <= 0.1 m (MaterialIDs 1) gets the
    rock properties (E, nu, Biot), i.e. rock up to R = 0.1 m; the mesh itself is unchanged
  - default case T0 = 300 degC only (spec v2); HM control (T0 = 15 degC) with --T0 15
  - full-field output at the snapshot times of spec Table 16: 2000, 2060, 3800, 5600 s

Original header (model v2):
Updated to the SAFENET-2 THM benchmark specification (Update-GREAT-THM.docx, Pan et al., 09.11.2026):
benchmark THM2 (intact disk, centre-point injection).

Changes with respect to version 1 (../THM_model/make_THM_prj.py):
  Table 2  matrix permeability k = 2.58e-19 m2 (Greywacke), 1e-19 m2 (Gneiss)   [v1: 1.409e-14]
           (E, nu, rho_s, porosity, Biot unchanged)
  Table 4  fluid compressibility c_f = 0                                         [v1: 4.4e-10 1/Pa]
           rho_f = 1000 kg/m3, mu = 1e-3 Pa s (constant for T0 <= 100 degC; variable for T0 = 300 degC)
  Table 6  solid lambda_s = 2.5 / 3.0 W/(m K), c_s = 850 / 800 J/(kg K),
           alpha_s = 8e-6 / 7e-6 1/K (Greywacke / Gneiss)                        [v1: bulk 2.8, 2700*790]
  Table 7  water lambda_f = 0.6, c_f = 4180 J/(kg K), beta_f = 0 (variable at 300 degC), T_in = 15 degC
  Table 8  load ramp 0-1800 s, hold 1800-2000 s, Q linear 0 -> Q0 over 2000-2120 s, Q0 until 5600 s;
           Q0 = 0.1 ml/min = 1.667e-9 m3/s, Q2D = Q0/B with B = 0.2 m           [v1: 1 l/min per m, 4000-4050 s]
           dt = 10 s (0-2000 s), 1 s (2000-2120 s), 5 s (2120-5600 s)
  Table 9  uniform T0 = 90 or 300 degC (T_ref = T0), adiabatic perimeter, p0 = 0 (gauge), p = 0 on the
           perimeter, 15 degC at the inlet only during inflow                    [v1: heating ramp 293->363 K at
                                                                                  the boundary, p0 = 1 MPa]
  Each THM run is paired with an isothermal HM control (T0 = T_in = 15 degC, same mesh/time steps/flow).

Case IDs: THM2-{GW|GN}-{load}-T{90|300}; HM control: HM2-{GW|GN}-{load}-T15

Usage: python make_THM2_scII_prj.py [--load A] [--T0 300 [15 ...]]  -> <case_id>.prj for both rocks
"""
import argparse

LOADS = {  # Table 3: PEE i & ia, i = 1..8 [MPa]
    "A": [10.0, 6.64, 4.46, 1.17, 1.0, 3.82, 7.80, 9.95],
    "E": [9.95, 10.0, 6.64, 4.46, 1.17, 1.0, 3.82, 7.80],
    "B": [7.80, 9.95, 10.0, 6.64, 4.46, 1.17, 1.0, 3.82],
    "F": [3.82, 7.80, 9.95, 10.0, 6.64, 4.46, 1.17, 1.0],
    "C": [1.0, 3.82, 7.80, 9.95, 10.0, 6.64, 4.46, 1.17],
}
# Table 2 + Table 6
ROCKS = {
    "GW": dict(name="greywacke", E=26.85e9, nu=0.27, rho=2650.0, phi=0.005, alpha=0.8, k=2.58e-19,
               lam_s=2.50, c_s=850.0, alpha_T=8.0e-6),
    "GN": dict(name="gneiss", E=83.9e9, nu=0.21, rho=2750.0, phi=0.001, alpha=0.6, k=1.0e-19,
               lam_s=3.00, c_s=800.0, alpha_T=7.0e-6),
}

# Table 4 / Table 7
WATER = dict(rho=1000.0, mu=1.0e-3, c=4180.0, lam=0.6)
T_IN_C = 15.0
B = 0.2                       # specimen thickness [m]
Q0 = 0.10e-6 / 60.0           # 0.1 ml/min [m3/s]
Q2D = Q0 / B                  # [m3/s per m]
T_RAMP_LOAD, T_INJ0, T_INJ1, T_END = 1800.0, 2000.0, 2120.0, 5600.0

# Variable water properties for T0 > 100 degC (Table 4/7: "variable"), compressed liquid at
# p_abs ~ p + 8.6 MPa (Table 1); fits to IAPWS-IF97 data, 15...300 degC, x = T - 273.15 [degC]
RHO_POLY = [-4.92318707e-08, 2.66626666e-05, -7.07790348e-03, 8.13185288e-02, 1.00362286e+03]  # max err 0.1 %
CP_POLY = [8.82625568e-07, -4.00672296e-04, 7.07378018e-02, -4.06898771e+00, 4.22046309e+03]   # max err 0.4 %
# viscosity: Vogel equation mu = 1e-3 exp(-3.7188 + 578.919/(T - 137.546)) [Pa s], T [K] (err < 7 %)


def poly(c, x="(temperature-273.15)"):
    n = len(c) - 1
    return " + ".join(f"({v:.10g})*{x}^{n - i}" if n - i > 0 else f"({v:.10g})" for i, v in enumerate(c))


def dpoly(c, x="(temperature-273.15)"):
    n = len(c) - 1
    return " + ".join(f"({(n - i) * v:.10g})*{x}^{n - i - 1}" if n - i > 1 else f"({v:.10g})"
                      for i, v in enumerate(c[:-1]))


def prop_const(name, v):
    return f"""                        <property>
                            <name>{name}</name>
                            <type>Constant</type>
                            <value>{v:g}</value>
                        </property>"""


def prop_func(name, expr, dexpr):
    return f"""                        <property>
                            <name>{name}</name>
                            <type>Function</type>
                            <value>
                                <expression>{expr}</expression>
                            </value>
                            <dvalue>
                                <variable_name>temperature</variable_name>
                                <expression>{dexpr}</expression>
                            </dvalue>
                            <dvalue>
                                <variable_name>liquid_phase_pressure</variable_name>
                                <expression>0</expression>
                            </dvalue>
                        </property>"""


def liquid(variable):
    w = WATER
    if not variable:
        return "\n".join([prop_const("density", w["rho"]), prop_const("viscosity", w["mu"]),
                          prop_const("specific_heat_capacity", w["c"]),
                          prop_const("thermal_conductivity", w["lam"]),
                          prop_const("thermal_expansivity", 0.0)])
    rho, drho = poly(RHO_POLY), dpoly(RHO_POLY)
    mu = "1e-3*exp(-3.7188 + 578.919/(temperature - 137.546))"
    dmu = f"-578.919/((temperature - 137.546)^2)*{mu}"
    return "\n".join([
        prop_func("density", rho, drho),
        prop_func("viscosity", mu, dmu),
        prop_func("specific_heat_capacity", poly(CP_POLY), dpoly(CP_POLY)),
        prop_const("thermal_conductivity", w["lam"]),
        # volumetric thermal expansion beta_f = -1/rho drho/dT (derivative w.r.t. T neglected)
        prop_func("thermal_expansivity", f"-({drho})/({rho})", "0"),
    ])


def medium(mid, r, biot, variable):
    return f"""        <medium id="{mid}">
            <phases>
                <phase>
                    <type>AqueousLiquid</type>
                    <properties>
{liquid(variable)}
                    </properties>
                </phase>
                <phase>
                    <type>Solid</type>
                    <properties>
{prop_const("density", r['rho'])}
{prop_const("thermal_conductivity", r['lam_s'])}
{prop_const("specific_heat_capacity", r['c_s'])}
{prop_const("thermal_expansivity", r['alpha_T'])}
                    </properties>
                </phase>
            </phases>
            <properties>
                <property>
                    <name>permeability</name>
                    <type>Constant</type>
                    <value>{r['k']:g}</value>
                </property>
                <property>
                    <name>porosity</name>
                    <type>Constant</type>
                    <value>{r['phi']:g}</value>
                </property>
                <property>
                    <name>biot_coefficient</name>
                    <type>Constant</type>
                    <value>{biot:g}</value>
                </property>
                <property>
                    <name>thermal_conductivity</name>
                    <type>EffectiveThermalConductivityPorosityMixing</type>
                </property>
            </properties>
        </medium>"""


def case_id(rock, load, T0c):
    return f"{'HM' if T0c == T_IN_C else 'THM'}2-{rock}-{load}-T{T0c:g}"


def prj(rock, load, T0c):
    r = ROCKS[rock]
    cid = case_id(rock, load, T0c)
    variable = T0c > 100.0
    T0, T_in = T0c + 273.15, T_IN_C + 273.15
    pee = LOADS[load]
    dss = [(pee[i] + pee[(i + 1) % 8]) / 2 for i in range(8)]
    meshes = ["domain"] + [f"PEE_{i}" for i in range(1, 9)] + [f"DSS_{i}" for i in range(1, 9)] + \
             ["outer_boundary", "fix_bottom", "fix_left", "center", "probe_circle"]
    mesh_xml = "\n".join(f"        <mesh>{m}.vtu</mesh>" for m in meshes)
    par, bc = [], []
    for kind, vals in (("PEE", pee), ("DSS", dss)):
        for i, v in enumerate(vals, 1):
            par.append(f"""        <parameter>
            <name>sigma_{kind}_{i}</name>
            <type>CurveScaled</type>
            <curve>load_ramp</curve>
            <parameter>sigma_{kind}_{i}_target</parameter>
        </parameter>
        <parameter>
            <name>sigma_{kind}_{i}_target</name>
            <type>Constant</type>
            <value>{v * 1e6:.6g}</value>
        </parameter>""")
            bc.append(f"""                <boundary_condition>
                    <mesh>{kind}_{i}</mesh>
                    <type>NormalTraction</type>
                    <parameter>sigma_{kind}_{i}</parameter>
                </boundary_condition>""")
    table = "\n".join(f"      PEE {i} & {i}a : {p:5.2f} MPa    DSS {i} & {i}a : {d:6.3f} MPa"
                      for i, (p, d) in enumerate(zip(pee, dss), 1))
    qm = Q2D * WATER["rho"]
    return f"""<?xml version="1.0" encoding="ISO-8859-1"?>
<!--
  GREAT cell THM benchmark, Scenario II of the GREAT Cell Modelling Specification v2 (22.11.2026)
  Case {cid}: {r['name']}, Load {load}, T0 = {T0c:g} degC, 2D plane strain (THERMO_HYDRO_MECHANICS)
  {'Isothermal HM control (T0 = T_in = 15 degC)' if T0c == T_IN_C else 'THM run, paired HM control: ' + case_id(rock, load, T_IN_C)}

  Mesh   : make_great_cell_mesh.py (quadratic triangles; T, p order 1, u order 2)
           MaterialIDs 0 (r <= 0.097 m) and 1 (0.097...0.100 m) both rock: no rubber sheath
  Schedule (Table 8):
     0 - 1800 s : linear ramp of Load {load} (Table 3, DSS = mean of the neighbouring PEE)
{table}
  1800 - 2000 s : hold, Q = 0
  2000 - 2120 s : Q linear 0 -> Q0, T_in = 15 degC at the centre node
  2120 - 5600 s : Q = Q0, T_in = 15 degC
     Q0 = 0.1 ml/min = {Q0:.4g} m3/s; Q2D = Q0/B = {Q2D:.4g} m3/s per m (B = {B} m)
     -> mass rate {qm:.4g} kg/s per m (equation_balance_type = mass)
  BC/IC (Table 9): T0 = {T0:.2f} K uniform (T_ref = T0), adiabatic perimeter;
     p0 = 0 (gauge), p = 0 on the perimeter; u = 0, zero initial stress
  Rock (Table 2/6): E = {r['E']/1e9:g} GPa, nu = {r['nu']}, Biot = {r['alpha']}, porosity = {r['phi']},
     rho_s = {r['rho']:g} kg/m3, k = {r['k']:g} m2, lambda_s = {r['lam_s']:g} W/(m K),
     c_s = {r['c_s']:g} J/(kg K), alpha_s = {r['alpha_T']:g} 1/K
  Water (Table 4/7): {'rho = 1000 kg/m3, mu = 1e-3 Pa s, c = 4180 J/(kg K), beta_f = 0, c_f = 0' if not variable else 'rho(T), mu(T), c(T) from IAPWS-IF97 fits (p_abs ~ 8.6 MPa), beta_f = -1/rho drho/dT, c_f = 0'},
     lambda_f = 0.6 W/(m K)
  No rubber sheath (spec v2, Scenario II): MaterialIDs 1 has the rock properties
-->
<OpenGeoSysProject>
    <meshes>
{mesh_xml}
    </meshes>
    <processes>
        <process>
            <name>THM</name>
            <type>THERMO_HYDRO_MECHANICS</type>
            <integration_order>4</integration_order>
            <constitutive_relation id="0">
                <type>LinearElasticIsotropic</type>
                <youngs_modulus>E_rock</youngs_modulus>
                <poissons_ratio>nu_rock</poissons_ratio>
            </constitutive_relation>
            <constitutive_relation id="1">
                <type>LinearElasticIsotropic</type>
                <youngs_modulus>E_rock</youngs_modulus>
                <poissons_ratio>nu_rock</poissons_ratio>
            </constitutive_relation>
            <process_variables>
                <temperature>temperature</temperature>
                <pressure>pressure</pressure>
                <displacement>displacement</displacement>
            </process_variables>
            <equation_balance_type>mass</equation_balance_type>
            <specific_body_force>0 0</specific_body_force>
            <secondary_variables>
                <secondary_variable internal_name="sigma" output_name="sigma"/>
                <secondary_variable internal_name="epsilon" output_name="epsilon"/>
                <secondary_variable internal_name="velocity" output_name="darcy_velocity"/>
            </secondary_variables>
        </process>
    </processes>
    <media>
{medium(0, r, r['alpha'], variable)}
{medium(1, r, r['alpha'], variable)}
    </media>
    <time_loop>
        <processes>
            <process ref="THM">
                <nonlinear_solver>basic_newton</nonlinear_solver>
                <convergence_criterion>
                    <type>PerComponentDeltaX</type>
                    <norm_type>NORM2</norm_type>
                    <abstols>1e-3 1e1 1e-8 1e-8</abstols>
                </convergence_criterion>
                <time_discretization>
                    <type>BackwardEuler</type>
                </time_discretization>
                <time_stepping>
                    <type>FixedTimeStepping</type>
                    <t_initial>0</t_initial>
                    <t_end>{T_END:g}</t_end>
                    <timesteps>
                        <!-- load ramp and hold 0 ... 2000 s -->
                        <pair>
                            <repeat>200</repeat>
                            <delta_t>10</delta_t>
                        </pair>
                        <!-- injection ramp 2000 ... 2120 s -->
                        <pair>
                            <repeat>120</repeat>
                            <delta_t>1</delta_t>
                        </pair>
                        <!-- constant injection 2120 ... 5600 s -->
                        <pair>
                            <repeat>696</repeat>
                            <delta_t>5</delta_t>
                        </pair>
                    </timesteps>
                </time_stepping>
            </process>
        </processes>
        <outputs>
            <output>
                <type>VTK</type>
                <prefix>{cid}_{{:meshname}}</prefix>
                <meshes>
                    <mesh>domain</mesh>
                </meshes>
                <timesteps>
                    <pair>
                        <repeat>1</repeat>
                        <each_steps>100000</each_steps>
                    </pair>
                </timesteps>
                <fixed_output_times>2000 2060 3800 5600</fixed_output_times>
                <variables>
                    <variable>temperature</variable>
                    <variable>pressure</variable>
                    <variable>displacement</variable>
                    <variable>sigma</variable>
                    <variable>epsilon</variable>
                    <variable>darcy_velocity</variable>
                </variables>
                <suffix>_ts_{{:timestep}}_t_{{:time}}</suffix>
            </output>
            <output>
                <type>VTK</type>
                <prefix>{cid}_{{:meshname}}</prefix>
                <meshes>
                    <mesh>center</mesh>
                    <mesh>probe_circle</mesh>
                </meshes>
                <timesteps>
                    <pair>
                        <repeat>1</repeat>
                        <each_steps>1</each_steps>
                    </pair>
                </timesteps>
                <variables>
                    <variable>temperature</variable>
                    <variable>pressure</variable>
                    <variable>displacement</variable>
                    <variable>sigma</variable>
                    <variable>epsilon</variable>
                </variables>
                <suffix>_ts_{{:timestep}}_t_{{:time}}</suffix>
            </output>
        </outputs>
    </time_loop>
    <parameters>
        <parameter>
            <name>E_rock</name>
            <type>Constant</type>
            <value>{r['E']:g}</value>
        </parameter>
        <parameter>
            <name>nu_rock</name>
            <type>Constant</type>
            <value>{r['nu']}</value>
        </parameter>
        <parameter>
            <name>T0</name>
            <type>Constant</type>
            <value>{T0:.2f}</value>
        </parameter>
        <parameter>
            <name>u0</name>
            <type>Constant</type>
            <values>0 0</values>
        </parameter>
        <parameter>
            <name>zero</name>
            <type>Constant</type>
            <value>0</value>
        </parameter>
        <parameter>
            <name>T_inj</name>
            <type>Constant</type>
            <value>{T_in:.2f}</value>
        </parameter>
        <parameter>
            <name>Q_inj</name>
            <type>Function</type>
            <!-- Table 8: Q2D = {Q2D:.4g} m3/s per m -> mass rate rho_w * Q2D; linear ramp 2000-2120 s -->
            <expression>(t &lt;= {T_INJ0:g}) ? 0 : ((t &lt; {T_INJ1:g}) ? {qm:.6g}*(t-{T_INJ0:g})/{T_INJ1 - T_INJ0:g} : {qm:.6g})</expression>
        </parameter>
{chr(10).join(par)}
    </parameters>
    <curves>
        <curve>
            <name>load_ramp</name>
            <coords>0 {T_RAMP_LOAD:g} {T_END:g}</coords>
            <values>0 1 1</values>
        </curve>
    </curves>
    <process_variables>
        <process_variable>
            <name>temperature</name>
            <components>1</components>
            <order>1</order>
            <initial_condition>T0</initial_condition>
            <!-- adiabatic outer perimeter (no BC); 15 degC at the inlet only during inflow -->
            <boundary_conditions>
                <boundary_condition>
                    <mesh>center</mesh>
                    <type>DirichletWithinTimeInterval</type>
                    <parameter>T_inj</parameter>
                    <time_interval>
                        <start>{T_INJ0 + 0.5:g}</start>
                        <end>{T_END:g}</end>
                    </time_interval>
                </boundary_condition>
            </boundary_conditions>
        </process_variable>
        <process_variable>
            <name>pressure</name>
            <components>1</components>
            <order>1</order>
            <initial_condition>zero</initial_condition>
            <boundary_conditions>
                <boundary_condition>
                    <mesh>outer_boundary</mesh>
                    <type>Dirichlet</type>
                    <parameter>zero</parameter>
                </boundary_condition>
            </boundary_conditions>
            <source_terms>
                <source_term>
                    <mesh>center</mesh>
                    <type>Nodal</type>
                    <parameter>Q_inj</parameter>
                </source_term>
            </source_terms>
        </process_variable>
        <process_variable>
            <name>displacement</name>
            <components>2</components>
            <order>2</order>
            <initial_condition>u0</initial_condition>
            <boundary_conditions>
                <boundary_condition>
                    <mesh>fix_bottom</mesh>
                    <type>Dirichlet</type>
                    <component>0</component>
                    <parameter>zero</parameter>
                </boundary_condition>
                <boundary_condition>
                    <mesh>fix_bottom</mesh>
                    <type>Dirichlet</type>
                    <component>1</component>
                    <parameter>zero</parameter>
                </boundary_condition>
                <boundary_condition>
                    <mesh>fix_left</mesh>
                    <type>Dirichlet</type>
                    <component>1</component>
                    <parameter>zero</parameter>
                </boundary_condition>
{chr(10).join(bc)}
            </boundary_conditions>
        </process_variable>
    </process_variables>
    <nonlinear_solvers>
        <nonlinear_solver>
            <name>basic_newton</name>
            <type>Newton</type>
            <max_iter>20</max_iter>
            <linear_solver>general_linear_solver</linear_solver>
        </nonlinear_solver>
    </nonlinear_solvers>
    <linear_solvers>
        <linear_solver>
            <name>general_linear_solver</name>
            <eigen>
                <solver_type>SparseLU</solver_type>
                <scaling>true</scaling>
            </eigen>
        </linear_solver>
    </linear_solvers>
</OpenGeoSysProject>
"""


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--load", default="A", choices=list(LOADS))
    ap.add_argument("--T0", type=float, nargs="+", default=[300.0], help="initial temperature(s) [degC]")
    a = ap.parse_args()
    for rock in ROCKS:
        for T0c in a.T0:
            fn = f"{case_id(rock, a.load, T0c)}.prj"
            open(fn, "w").write(prj(rock, a.load, T0c))
            print("wrote", fn)
