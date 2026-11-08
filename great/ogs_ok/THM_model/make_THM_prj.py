"""
OGS THM model (THERMO_HYDRO_MECHANICS, 2D plane strain) of the GREAT cell sample, Load case A.
It combines:
  - the M model (benchmark M1, Mollaali et al. 2026, JRMGE, doi:10.1016/j.jrmge.2026.03.060):
    geometry, rubber sheath, PEE/DSS loads (Table 4, Eq. 19), elasticity and Biot/porosity (Table 1),
    fluid compressibility (Table 3), rigid-body constraints (Eq. 18)
  - the TH model (TH_granite_injection.prj, STIMTEC): heating ramp on the outer boundary,
    water injection at the centre 4000-4050 s, shut-in until 5000 s, permeability, p0, thermal
    properties (bulk values identical to the TH/T models)

Usage: python make_THM_prj.py [--load A] [--alpha_T 8e-6]  -> THM_<rock>_load<X>.prj
"""
import argparse

LOADS = {  # Table 4: PEE i & ia, i = 1..8 [MPa]
    "A": [10.0, 6.64, 4.46, 1.17, 1.0, 3.82, 7.80, 9.95],
    "E": [9.95, 10.0, 6.64, 4.46, 1.17, 1.0, 3.82, 7.80],
    "B": [7.80, 9.95, 10.0, 6.64, 4.46, 1.17, 1.0, 3.82],
    "F": [3.82, 7.80, 9.95, 10.0, 6.64, 4.46, 1.17, 1.0],
    "C": [1.0, 3.82, 7.80, 9.95, 10.0, 6.64, 4.46, 1.17],
}
# Table 1: E [Pa], nu, rho_s [kg/m3], porosity, Biot coefficient
ROCKS = {
    "greywacke": dict(E=26.85e9, nu=0.27, rho=2650.0, phi=0.005, alpha=0.8),
    "gneiss": dict(E=83.9e9, nu=0.21, rho=2750.0, phi=0.001, alpha=0.6),
}
SHEATH = dict(E=0.1e9, nu=0.4, rho=1200.0, alpha=0.0)   # Table 1 (Biot 0)

# from the TH model
K_PERM = 1.409e-14          # m2 (STIMTEC calibration V1)
P0, T0 = 1.0e6, 293.0
Q_INJ = 1.6666667e-5        # m3/s per m thickness (1 l/min)
LAMBDA_BULK, RHOC_BULK = 2.8, 2700 * 790.0
WATER = dict(rho=1000.0, mu=1.0e-3, c=4180.0, lam=0.6, beta=4.4e-10)  # beta: Table 3


def solid_thermal(phi, rho_s):
    """solid lambda and c such that the bulk values equal those of the TH model"""
    lam_s = (LAMBDA_BULK - phi * WATER["lam"]) / (1 - phi)
    c_s = (RHOC_BULK - phi * WATER["rho"] * WATER["c"]) / ((1 - phi) * rho_s)
    return lam_s, c_s


def medium(mid, phi, alpha, rho_s, alpha_T):
    lam_s, c_s = solid_thermal(phi, rho_s)
    w = WATER
    return f"""        <medium id="{mid}">
            <phases>
                <phase>
                    <type>AqueousLiquid</type>
                    <properties>
                        <property>
                            <name>density</name>
                            <type>Linear</type>
                            <reference_value>{w['rho']:g}</reference_value>
                            <independent_variable>
                                <variable_name>liquid_phase_pressure</variable_name>
                                <reference_condition>{P0:g}</reference_condition>
                                <slope>{w['beta']:g}</slope>
                            </independent_variable>
                        </property>
                        <property>
                            <name>viscosity</name>
                            <type>Constant</type>
                            <value>{w['mu']:g}</value>
                        </property>
                        <property>
                            <name>specific_heat_capacity</name>
                            <type>Constant</type>
                            <value>{w['c']:g}</value>
                        </property>
                        <property>
                            <name>thermal_conductivity</name>
                            <type>Constant</type>
                            <value>{w['lam']:g}</value>
                        </property>
                        <property>
                            <name>thermal_expansivity</name>
                            <type>Constant</type>
                            <value>0</value>
                        </property>
                    </properties>
                </phase>
                <phase>
                    <type>Solid</type>
                    <properties>
                        <property>
                            <name>density</name>
                            <type>Constant</type>
                            <value>{rho_s:g}</value>
                        </property>
                        <property>
                            <name>thermal_conductivity</name>
                            <type>Constant</type>
                            <value>{lam_s:.6g}</value>
                        </property>
                        <property>
                            <name>specific_heat_capacity</name>
                            <type>Constant</type>
                            <value>{c_s:.6g}</value>
                        </property>
                        <property>
                            <name>thermal_expansivity</name>
                            <type>Constant</type>
                            <value>{alpha_T:g}</value>
                        </property>
                    </properties>
                </phase>
            </phases>
            <properties>
                <property>
                    <name>permeability</name>
                    <type>Constant</type>
                    <value>{K_PERM:g}</value>
                </property>
                <property>
                    <name>porosity</name>
                    <type>Constant</type>
                    <value>{phi:g}</value>
                </property>
                <property>
                    <name>biot_coefficient</name>
                    <type>Constant</type>
                    <value>{alpha:g}</value>
                </property>
                <property>
                    <name>thermal_conductivity</name>
                    <type>EffectiveThermalConductivityPorosityMixing</type>
                </property>
            </properties>
        </medium>"""


def prj(rock, load, alpha_T):
    r = ROCKS[rock]
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
            <type>Constant</type>
            <value>{v * 1e6:.6g}</value>
        </parameter>""")
            bc.append(f"""                <boundary_condition>
                    <mesh>{kind}_{i}</mesh>
                    <type>NormalTraction</type>
                    <parameter>sigma_{kind}_{i}</parameter>
                </boundary_condition>""")
    lam_s, c_s = solid_thermal(r["phi"], r["rho"])
    lam_sh, c_sh = solid_thermal(r["phi"], SHEATH["rho"])
    table = "\n".join(f"      PEE {i} & {i}a : {p:5.2f} MPa    DSS {i} & {i}a : {d:6.3f} MPa"
                      for i, (p, d) in enumerate(zip(pee, dss), 1))
    s0 = r["alpha"] * P0
    return f"""<?xml version="1.0" encoding="ISO-8859-1"?>
<!--
  GREAT cell THM model - {rock}, Load {load}, 2D plane strain (THERMO_HYDRO_MECHANICS)
  = M model (benchmark M1, Mollaali et al. 2026, JRMGE, doi:10.1016/j.jrmge.2026.03.060)
  + TH model (TH_granite_injection.prj: heating, injection, shut-in)

  Mesh   : make_great_cell_mesh.py (quadratic triangles; T, p order 1, u order 2;
           h = 1.5 mm at the centre -> 6 mm at r = 0.065 m, 2 mm in the sheath)
           rock r <= 0.097 m (MaterialIDs 0), rubber sheath 0.097...0.100 m (MaterialIDs 1)
  Phase 0 (t = 0)          : mechanical Load {load} applied in the first time step
                             (normal stress on PEE's, Table 4; DSS's = mean of neighbours, Eq. 19):
{table}
  Phase 1 (0 ... 4000 s)   : outer boundary T ramp 293 -> 363 K (0...325 s), then 363 K
  Phase 2 (4000 ... 4050 s): injection at the centre node, Q = 1 l/min = 1.6667e-5 m3/s per m
                             (mass rate 1.6667e-2 kg/s per m, equation_balance_type = mass),
                             T_inj = 293 K (Dirichlet on the centre node)
  Phase 3 (4050 ... 5000 s): shut-in
  Hydraulic BC : p = p0 = 1 MPa on the outer boundary (drained)
  Initial      : T0 = 293 K, p0 = 1 MPa, u = 0, initial effective stress alpha*p0 so that the
                 initial state (p = p0, no load) is stress free in terms of total stress

  Rock ({rock}, Table 1): E = {r['E']/1e9:g} GPa, nu = {r['nu']}, Biot alpha = {r['alpha']},
                 porosity = {r['phi']}, rho_s = {r['rho']:g} kg/m3
  Sheath (Table 1): E = 0.1 GPa, nu = 0.4, Biot alpha = 0
  Fluid  : rho = 1000 kg/m3 at p0, compressibility 4.4e-10 1/Pa (Table 3), mu = 1e-3 Pa s,
           c = 4180 J/(kg K), lambda = 0.6 W/(m K), no thermal expansion (as in the TH model)
  From TH model: k = 1.409e-14 m2 (rock and sheath); bulk lambda = 2.8 W/(m K) and
           bulk rho*c = 2700*790 J/(m3 K) in rock and sheath
           -> rock  : lambda_s = {lam_s:.5g}, c_s = {c_s:.5g}
           -> sheath: lambda_s = {lam_sh:.5g}, c_s = {c_sh:.5g}
  Assumptions: linear thermal expansion coefficient of the solid {alpha_T:g} 1/K (rock and sheath);
           sheath hydraulically and thermally like the rock (drained boundary as in the TH model);
           PEE arc 18 deg, DSS arc 4.5 deg; no gravity
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
                <youngs_modulus>E_sheath</youngs_modulus>
                <poissons_ratio>nu_sheath</poissons_ratio>
            </constitutive_relation>
            <process_variables>
                <temperature>temperature</temperature>
                <pressure>pressure</pressure>
                <displacement>displacement</displacement>
            </process_variables>
            <!-- mass balance: needed for the pressure-dependent fluid density (compressibility);
                 the injection source term is therefore a mass rate [kg/s per m] -->
            <equation_balance_type>mass</equation_balance_type>
            <initial_stress>sigma0</initial_stress>
            <specific_body_force>0 0</specific_body_force>
            <secondary_variables>
                <secondary_variable internal_name="sigma" output_name="sigma"/>
                <secondary_variable internal_name="epsilon" output_name="epsilon"/>
                <secondary_variable internal_name="velocity" output_name="darcy_velocity"/>
            </secondary_variables>
        </process>
    </processes>
    <media>
{medium(0, r['phi'], r['alpha'], r['rho'], alpha_T)}
{medium(1, r['phi'], SHEATH['alpha'], SHEATH['rho'], alpha_T)}
    </media>
    <time_loop>
        <processes>
            <process ref="THM">
                <nonlinear_solver>basic_newton</nonlinear_solver>
                <convergence_criterion>
                    <type>PerComponentDeltaX</type>
                    <norm_type>NORM2</norm_type>
                    <abstols>1e-1 1e1 1e-7 1e-7</abstols>
                </convergence_criterion>
                <time_discretization>
                    <type>BackwardEuler</type>
                </time_discretization>
                <time_stepping>
                    <type>FixedTimeStepping</type>
                    <t_initial>0</t_initial>
                    <t_end>5000</t_end>
                    <timesteps>
                        <!-- ramp of the boundary temperature 0 ... 400 s -->
                        <pair>
                            <repeat>80</repeat>
                            <delta_t>5</delta_t>
                        </pair>
                        <!-- heating 400 ... 4000 s -->
                        <pair>
                            <repeat>90</repeat>
                            <delta_t>40</delta_t>
                        </pair>
                        <!-- injection 4000 ... 4050 s -->
                        <pair>
                            <repeat>50</repeat>
                            <delta_t>1</delta_t>
                        </pair>
                        <!-- early shut-in 4050 ... 4150 s -->
                        <pair>
                            <repeat>50</repeat>
                            <delta_t>2</delta_t>
                        </pair>
                        <!-- recovery 4150 ... 5000 s -->
                        <pair>
                            <repeat>85</repeat>
                            <delta_t>10</delta_t>
                        </pair>
                    </timesteps>
                </time_stepping>
            </process>
        </processes>
        <outputs>
            <output>
                <type>VTK</type>
                <prefix>THM_{rock}_load{load}_{{:meshname}}</prefix>
                <meshes>
                    <mesh>domain</mesh>
                </meshes>
                <timesteps>
                    <pair>
                        <repeat>1</repeat>
                        <each_steps>100000</each_steps>
                    </pair>
                </timesteps>
                <fixed_output_times>5 325 2000 4000 4001 4005 4010 4025 4050 4060 4100 4200 4500 5000</fixed_output_times>
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
                <prefix>THM_{rock}_load{load}_{{:meshname}}</prefix>
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
            <name>E_sheath</name>
            <type>Constant</type>
            <value>{SHEATH['E']:g}</value>
        </parameter>
        <parameter>
            <name>nu_sheath</name>
            <type>Constant</type>
            <value>{SHEATH['nu']}</value>
        </parameter>
        <!-- initial effective stress alpha*p0 (rock) / 0 (sheath, alpha = 0): zero total stress at p = p0 -->
        <parameter>
            <name>sigma0</name>
            <type>Group</type>
            <group_id_property>MaterialIDs</group_id_property>
            <index_values>
                <index>0</index>
                <values>{s0:g} {s0:g} {s0:g} 0</values>
            </index_values>
            <index_values>
                <index>1</index>
                <values>0 0 0 0</values>
            </index_values>
        </parameter>
        <parameter>
            <name>T0</name>
            <type>Constant</type>
            <value>{T0:g}</value>
        </parameter>
        <parameter>
            <name>p0</name>
            <type>Constant</type>
            <value>{P0:g}</value>
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
            <name>one</name>
            <type>Constant</type>
            <value>1</value>
        </parameter>
        <parameter>
            <name>T_outer</name>
            <type>CurveScaled</type>
            <curve>T_outer_curve</curve>
            <parameter>one</parameter>
        </parameter>
        <parameter>
            <name>T_inj</name>
            <type>Constant</type>
            <value>293</value>
        </parameter>
        <parameter>
            <name>Q_inj</name>
            <type>Function</type>
            <!-- 1 l/min = {Q_INJ:g} m3/s per m  ->  mass rate rho_w * Q -->
            <expression>(t &gt; 4000 and t &lt;= 4050) ? {Q_INJ * WATER['rho']:g} : 0</expression>
        </parameter>
{chr(10).join(par)}
    </parameters>
    <curves>
        <curve>
            <name>T_outer_curve</name>
            <coords>0 325 5000</coords>
            <values>293 363 363</values>
        </curve>
    </curves>
    <process_variables>
        <process_variable>
            <name>temperature</name>
            <components>1</components>
            <order>1</order>
            <initial_condition>T0</initial_condition>
            <boundary_conditions>
                <boundary_condition>
                    <mesh>outer_boundary</mesh>
                    <type>Dirichlet</type>
                    <parameter>T_outer</parameter>
                </boundary_condition>
                <boundary_condition>
                    <mesh>center</mesh>
                    <type>DirichletWithinTimeInterval</type>
                    <parameter>T_inj</parameter>
                    <time_interval>
                        <start>4000.5</start>
                        <end>4050</end>
                    </time_interval>
                </boundary_condition>
            </boundary_conditions>
        </process_variable>
        <process_variable>
            <name>pressure</name>
            <components>1</components>
            <order>1</order>
            <initial_condition>p0</initial_condition>
            <boundary_conditions>
                <boundary_condition>
                    <mesh>outer_boundary</mesh>
                    <type>Dirichlet</type>
                    <parameter>p0</parameter>
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
    ap.add_argument("--alpha_T", type=float, default=8e-6, help="solid linear thermal expansion [1/K]")
    a = ap.parse_args()
    for rock in ROCKS:
        fn = f"THM_{rock}_load{a.load}.prj"
        open(fn, "w").write(prj(rock, a.load, a.alpha_T))
        print("wrote", fn)
