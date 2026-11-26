"""Export of the OGS THM2 results, Scenario II, in the format of the
GREAT Cell Modelling Specification v2 (22.11.2026), Section 8.2 (Tables 13-16).

Per case <case> (e.g. THM2-GW-A-T300) it writes to <results_dir>:
  UFZ_<case>_points.csv        four-point histories P0, P90, P180, P270 on r_p = 0.065 m, every time step
                               t, T_P0, p_P0, eps_vol_P0, ..., T_P270, p_P270, eps_vol_P270
  UFZ_<case>_inlet.csv         injection data at the centre node, every time step
                               t, p_in, T_in, Q_ml_min, V_cum_ml
  UFZ_<case>_profiles.csv      circumferential profiles on r_p at 2000, 2060, 3800, 5600 s,
                               theta = 0...359 deg at 1 deg spacing: t, theta_deg, T, p, eps_vol
  UFZ_<case>_field_t<time>.vtu full-field snapshots at 2000, 2060, 3800, 5600 s with
                               u_magnitude [m], sigma_trace [MPa] (total stress, 3D trace incl. sigma_zz),
                               sigma_eff_trace [MPa] (Biot effective stress), eps_vol [-], p [MPa],
                               q [m/s] (Darcy flux), T [degC], u [m]

Units (spec Sect. 8.2): t in s, T in degC, p in MPa (gauge), strain as fraction, Q in ml/min,
V in ml, all other quantities SI. Sign convention: tension positive.

Interpolation rule (identical for T, p and eps_vol): values at the corner nodes of the probe-circle
elements (every 5 deg; p and T are linear fields and are output at corner nodes only) are
interpolated linearly in theta. P0, P90, P180 and P270 coincide with corner nodes, so the
four-point histories are nodal values.

Usage: python export_THM2_scII.py <ogs_out_dir> <results_dir> [--cases THM2-GW-A-T300 THM2-GN-A-T300]
"""
import argparse
import glob
import os
import re
import numpy as np
import pyvista as pv

TK = 273.15
B, Q0 = 0.2, 0.1e-6 / 60.0            # specimen thickness [m], Q0 = 0.1 ml/min [m3/s]
SNAP = [2000.0, 2060.0, 3800.0, 5600.0]
BIOT = {"GW": 0.8, "GN": 0.6}         # Biot coefficient (spec Table 4)
POINTS = {"P0": 0.0, "P90": 90.0, "P180": 180.0, "P270": 270.0}


def series(prefix):
    fs = glob.glob(f"{prefix}_ts_*_t_*.vtu")
    t = np.array([float(re.search(r"_t_([0-9.]+)\.vtu$", f).group(1)) for f in fs])
    o = np.argsort(t)
    return t[o], [fs[i] for i in o]


def q_of_t(t):
    """3D injection rate Q(t) [m3/s] of the schedule (spec Table 11)."""
    return np.where(t <= 2000, 0.0, np.where(t < 2120, Q0 * (t - 2000) / 120, Q0))


class Circle:
    """Probe circle: linear interpolation in theta between the element corner nodes (every 5 deg)."""

    def __init__(self, mesh):
        th = np.degrees(np.arctan2(mesh.points[:, 1], mesh.points[:, 0])) % 360
        cn = np.unique(mesh.cells_dict[21][:, :2])      # corner nodes of the quadratic edges
        o = np.argsort(th[cn])
        self.nodes, self.th = cn[o], th[cn][o]

    def interp(self, values, theta):
        v = values[self.nodes]
        return np.interp(np.asarray(theta, float) % 360, np.r_[self.th, self.th[0] + 360],
                         np.r_[v, v[0]])


def fill_midnodes(m, v):
    """Order-1 fields (p, T) are written at the corner nodes only; set the mid nodes of the
    quadratic triangles to the mean of their edge end nodes (linear interpolation)."""
    v = np.array(v, float)
    c = m.cells_dict[22]
    for mid, (a, b) in zip((3, 4, 5), ((0, 1), (1, 2), (2, 0))):
        v[c[:, mid]] = 0.5 * (v[c[:, a]] + v[c[:, b]])
    return v


def export(case, out, res):
    t, fs = series(f"{out}/{case}_probe_circle")
    if len(t) < 2:
        print(f"{case}: no results in {out}")
        return
    circ = Circle(pv.read(fs[0]))
    ang = np.array(list(POINTS.values()))
    prof_th = np.arange(0.0, 360.0, 1.0)
    P_rows, prof = [], []
    for ti, f in zip(t, fs):
        d = pv.read(f).point_data
        T = d["temperature"] - TK
        p = d["pressure"] / 1e6
        e = d["epsilon"][:, 0] + d["epsilon"][:, 1]           # plane strain: eps_zz = 0
        row = [ti]
        for v in zip(circ.interp(T, ang), circ.interp(p, ang), circ.interp(e, ang)):
            row += list(v)
        P_rows.append(row)
        if np.any(np.isclose(ti, SNAP)):
            for th, a, b, c in zip(prof_th, circ.interp(T, prof_th), circ.interp(p, prof_th),
                                   circ.interp(e, prof_th)):
                prof.append([ti, th, a, b, c])
    hdr = "t," + ",".join(f"T_{n},p_{n},eps_vol_{n}" for n in POINTS)
    np.savetxt(f"{res}/UFZ_{case}_points.csv", np.array(P_rows), delimiter=",", header=hdr,
               comments="", fmt="%.10g")
    np.savetxt(f"{res}/UFZ_{case}_profiles.csv", np.array(prof), delimiter=",",
               header="t,theta_deg,T,p,eps_vol", comments="", fmt="%.10g")

    # inlet: centre node
    tc, fc = series(f"{out}/{case}_center")
    Tin = np.array([pv.read(f).point_data["temperature"][0] for f in fc]) - TK
    pin = np.array([pv.read(f).point_data["pressure"][0] for f in fc]) / 1e6
    q = q_of_t(tc)
    V = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(tc))]) * 1e6   # ml
    np.savetxt(f"{res}/UFZ_{case}_inlet.csv", np.column_stack([tc, pin, Tin, q * 6e7, V]),
               delimiter=",", header="t,p_in,T_in,Q_ml_min,V_cum_ml", comments="", fmt="%.10g")

    # full-field snapshots
    td, fd = series(f"{out}/{case}_domain")
    for ts in SNAP:
        k = int(np.argmin(np.abs(td - ts)))
        if abs(td[k] - ts) > 1e-6:
            print(f"  warning: no domain output at t = {ts:g} s")
            continue
        m = pv.read(fd[k])
        g = pv.UnstructuredGrid(m.cells, m.celltypes, m.points)
        pd, cd = m.point_data, m.cell_data
        g.point_data["u_magnitude"] = np.linalg.norm(pd["displacement"], axis=1)
        # OGS THERMO_HYDRO_MECHANICS writes the effective stress sigma' (Kelvin vector xx, yy, zz, xy);
        # total stress sigma = sigma' - alpha_B p I (tension positive)
        sig = pd["sigma"]
        tr_eff = (sig[:, 0] + sig[:, 1] + sig[:, 2]) / 1e6
        p_full = fill_midnodes(m, pd["pressure"]) / 1e6
        g.point_data["sigma_trace"] = tr_eff - 3 * BIOT[case.split("-")[1]] * p_full
        g.point_data["sigma_eff_trace"] = tr_eff
        g.point_data["eps_vol"] = pd["epsilon"][:, 0] + pd["epsilon"][:, 1]
        g.point_data["p"] = fill_midnodes(m, pd["pressure"]) / 1e6
        g.point_data["T"] = fill_midnodes(m, pd["temperature"]) - TK
        g.point_data["u"] = pd["displacement"]
        if "darcy_velocity" in pd:
            g.point_data["q"] = pd["darcy_velocity"]
        elif "darcy_velocity" in cd:
            g.cell_data["q"] = cd["darcy_velocity"]
        if "MaterialIDs" in cd:
            g.cell_data["MaterialIDs"] = cd["MaterialIDs"]
        g.field_data["time_s"] = np.array([td[k]])
        g.save(f"{res}/UFZ_{case}_field_t{ts:g}.vtu")
    print(f"{case}: {len(t)} steps, t_end = {t[-1]:g} s, V_cum = {V[-1]:.3f} ml -> {res}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("res")
    ap.add_argument("--cases", nargs="+", default=["THM2-GW-A-T300", "THM2-GN-A-T300"])
    a = ap.parse_args()
    os.makedirs(a.res, exist_ok=True)
    for c in a.cases:
        export(c, a.out, a.res)
