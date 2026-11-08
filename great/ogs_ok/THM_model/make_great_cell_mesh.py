"""
2D plane-strain mesh of the GREAT cell sample for the OGS THM model (quadratic elements,
displacement order 2), derived from the M1 mesh script (benchmark M1,
Mollaali et al. 2026, JRMGE, doi:10.1016/j.jrmge.2026.03.060).

Geometry (paper, Sect. 3.1 / Fig. 4b)
  - rock disk           r <= R_ROCK  = 0.097 m   MaterialIDs = 0
  - rubber sheath       R_ROCK < r <= R_OUT = 0.100 m   MaterialIDs = 1
  - probe circle        r = R_PROBE = 0.065 m (embedded, volumetric strain vs angle)
  - centre node         (0, 0, 0)
  - outer boundary split into 16 PEE arcs (Pressure-Exerting Elements) and
    16 DSS arcs (Dynamic-Sealing Strips). PEE1 is centred at 0 deg, PEE i at
    -(i-1)*22.5 deg (clockwise), PEE ia opposite (+180 deg); DSS i lies between
    PEE i and PEE i+1 (DSS 8 between PEE 8 and PEE 1a).
    ASSUMPTION: PEE arc width W_PEE = 18 deg, DSS arc width 22.5 - W_PEE = 4.5 deg
    (not stated in the paper; estimated from Fig. 4b). Change with --w_pee.
  - rigid-body points (paper Eq. 18): P_bottom = (0, -0.09894) -> ux = uy = 0,
                                      P_left   = (-0.09894, 0) -> uy = 0

Output (VTU, OGS-ready, all subdomains with bulk_node_ids / bulk_element_ids)
  domain.vtu, PEE_1 ... PEE_8.vtu (each = arcs i and ia), DSS_1 ... DSS_8.vtu,
  outer_boundary.vtu, probe_circle.vtu, rock_disk.vtu, center.vtu,
  fix_bottom.vtu, fix_left.vtu

Usage: python make_great_cell_mesh.py [--lc 0.006] [--lc_sheath 0.002] [--lc_center 0.0015]
       [--order 2] [--w_pee 18] [--out .]
"""
import argparse
import os

import gmsh
import meshio
import numpy as np

R_OUT, R_ROCK, R_PROBE = 0.100, 0.097, 0.065
R_FIX = 0.09894
DTH = 22.5


def segment_angles(w_pee):
    """Return list of (name, start_deg, end_deg) for the 32 boundary arcs (ccw)."""
    hw = 0.5 * w_pee
    arcs = []
    for k in range(16):                       # k = 0..15 counter-clockwise from PEE1
        i = (-k) % 16                         # clockwise index from PEE1
        c = k * DTH                           # ccw centre angle of this PEE
        pee = i % 8 + 1
        name = f"PEE_{pee}"
        arcs.append((name, c - hw, c + hw))
        # DSS ccw-next of this PEE lies between PEE(i) and PEE(i-1) in clockwise numbering
        j = (i - 1) % 16
        dss = j % 8 + 1                       # DSS j sits between PEE j and PEE j+1
        arcs.append((f"DSS_{dss}", c + hw, c + DTH - hw))
    return arcs


def build(lc, lc_sheath, lc_center, w_pee, order, outdir):
    arcs = segment_angles(w_pee)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("great_cell_M1")
    g = gmsh.model.geo
    pc = g.addPoint(0, 0, 0, lc)

    def ring(R, angles, size):
        pts = [g.addPoint(R * np.cos(np.radians(a)), R * np.sin(np.radians(a)), 0, size)
               for a in angles]
        crv = [g.addCircleArc(pts[n], pc, pts[(n + 1) % len(pts)]) for n in range(len(pts))]
        return crv

    # outer boundary: break points at all PEE/DSS ends
    outer_angles = [a[1] for a in arcs]
    c_out = ring(R_OUT, outer_angles, lc_sheath)
    c_rock = ring(R_ROCK, [0, 90, 180, 270], lc_sheath)
    c_probe = ring(R_PROBE, [0, 90, 180, 270], lc)

    l_out = g.addCurveLoop(c_out)
    l_rock = g.addCurveLoop(c_rock)
    l_probe = g.addCurveLoop(c_probe)
    s_sheath = g.addPlaneSurface([l_out, l_rock])
    s_rock_out = g.addPlaneSurface([l_rock, l_probe])
    s_rock_in = g.addPlaneSurface([l_probe])
    p_bot = g.addPoint(0, -R_FIX, 0, lc_sheath)
    p_left = g.addPoint(-R_FIX, 0, 0, lc_sheath)
    g.synchronize()
    gmsh.model.mesh.embed(0, [pc], 2, s_rock_in)
    gmsh.model.mesh.embed(0, [p_bot, p_left], 2, s_sheath)

    # size: lc in rock, lc_sheath in the sheath and near the rock/sheath interface
    f = gmsh.model.mesh.field
    fd = f.add("Distance"); f.setNumbers(fd, "CurvesList", c_rock)
    ft = f.add("Threshold"); f.setNumber(ft, "InField", fd)
    f.setNumber(ft, "SizeMin", lc_sheath); f.setNumber(ft, "SizeMax", lc)
    f.setNumber(ft, "DistMin", R_OUT - R_ROCK); f.setNumber(ft, "DistMax", 0.01)
    # refinement towards the injection point (cooling front / pressure singularity)
    fc = f.add("Distance"); f.setNumbers(fc, "PointsList", [pc])
    ftc = f.add("Threshold"); f.setNumber(ftc, "InField", fc)
    f.setNumber(ftc, "SizeMin", lc_center); f.setNumber(ftc, "SizeMax", lc)
    f.setNumber(ftc, "DistMin", 0.0); f.setNumber(ftc, "DistMax", R_PROBE)
    fm = f.add("Min"); f.setNumbers(fm, "FieldsList", [ft, ftc])
    f.setAsBackgroundMesh(fm)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    gmsh.model.addPhysicalGroup(2, [s_rock_in, s_rock_out], 1, "rock")
    gmsh.model.addPhysicalGroup(2, [s_sheath], 2, "sheath")
    names = sorted({a[0] for a in arcs})
    for n, name in enumerate(names):
        gmsh.model.addPhysicalGroup(
            1, [c_out[k] for k, a in enumerate(arcs) if a[0] == name], 10 + n, name)
    gmsh.model.addPhysicalGroup(1, c_probe, 3, "probe_circle")
    gmsh.model.addPhysicalGroup(0, [pc], 4, "center")
    gmsh.model.addPhysicalGroup(0, [p_bot], 5, "fix_bottom")
    gmsh.model.addPhysicalGroup(0, [p_left], 6, "fix_left")
    gmsh.model.addPhysicalGroup(2, [s_rock_in], 7, "probe_disk")

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Laplace2D")
    if order == 2:
        gmsh.model.mesh.setOrder(2)
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()
    gmsh.write(os.path.join(outdir, "great_cell.msh"))

    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    o = np.argsort(tags); tags, xyz = tags[o], xyz[o]
    xyz[:, 2] = 0.0
    t2i = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
    t2i[tags.astype(np.int64)] = np.arange(len(tags))

    def elems(dim, phys):
        conn = []
        for ent in gmsh.model.getEntitiesForPhysicalGroup(dim, phys):
            types, _, ntags = gmsh.model.mesh.getElements(dim, ent)
            for et, nt in zip(types, ntags):
                nn = gmsh.model.mesh.getElementProperties(et)[3]
                conn.append(t2i[nt.astype(np.int64)].reshape(-1, nn))
        return np.vstack(conn)

    def node(phys):
        return int(t2i[gmsh.model.mesh.getNodesForPhysicalGroup(0, phys)[0][0]])

    tri_rock, tri_sheath = elems(2, 1), elems(2, 2)
    tri_probe_disk = elems(2, 7)
    segs = {name: elems(1, 10 + n) for n, name in enumerate(names)}
    seg_probe = elems(1, 3)
    n_c, n_b, n_l = node(4), node(5), node(6)
    gmsh.finalize()

    tri = np.vstack([tri_rock, tri_sheath])
    mat = np.r_[np.zeros(len(tri_rock)), np.ones(len(tri_sheath))].astype(np.int32)
    p = xyz[tri]
    area2 = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - \
            (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])
    flip = [0, 2, 1] if tri.shape[1] == 3 else [0, 2, 1, 5, 4, 3]
    tri[area2 < 0] = tri[area2 < 0][:, flip]
    tri_t = "triangle" if tri.shape[1] == 3 else "triangle6"
    line_t = "line" if tri.shape[1] == 3 else "line3"
    meshio.write(os.path.join(outdir, "domain.vtu"),
                 meshio.Mesh(xyz, [(tri_t, tri)], cell_data={"MaterialIDs": [mat]}),
                 binary=True)

    edge2tri = {}
    for e, t in enumerate(tri):
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):  # corner nodes
            edge2tri.setdefault(frozenset((a, b)), e)
    tri_key = {frozenset(t[:3]): e for e, t in enumerate(tri)}

    def write_sub(name, cells, celltype, bulk_eids=None):
        used = np.unique(cells)
        loc = {gid: i for i, gid in enumerate(used)}
        lc_ = np.vectorize(loc.get)(cells)
        cd = {} if bulk_eids is None else {"bulk_element_ids": [np.asarray(bulk_eids, np.uint64)]}
        m = meshio.Mesh(xyz[used], [(celltype, lc_)],
                        point_data={"bulk_node_ids": used.astype(np.uint64)}, cell_data=cd)
        meshio.write(os.path.join(outdir, name), m, binary=True)

    def seg_ids(s):
        return [edge2tri[frozenset(x[:2])] for x in s]

    for name, s in segs.items():
        write_sub(f"{name}.vtu", s, line_t, seg_ids(s))
    allseg = np.vstack(list(segs.values()))
    write_sub("outer_boundary.vtu", allseg, line_t, seg_ids(allseg))
    write_sub("probe_circle.vtu", seg_probe, line_t, seg_ids(seg_probe))
    write_sub("probe_disk.vtu", tri_probe_disk, tri_t,
              [tri_key[frozenset(t[:3])] for t in tri_probe_disk])
    for fname, n in (("center.vtu", n_c), ("fix_bottom.vtu", n_b), ("fix_left.vtu", n_l)):
        write_sub(fname, np.array([[n]]), "vertex")

    A = 0.5 * np.abs(area2)
    r = np.hypot(xyz[:, 0], xyz[:, 1])
    print(f"nodes {len(xyz)}, triangles {len(tri)} (rock {len(tri_rock)}, sheath {len(tri_sheath)})")
    print(f"area {A.sum():.6e} vs {np.pi*R_OUT**2:.6e}; rock {A[mat==0].sum():.6e} vs {np.pi*R_ROCK**2:.6e}")
    print(f"nodes on probe circle {np.sum(np.abs(r-R_PROBE)<1e-9)}, fix points {xyz[n_b]}, {xyz[n_l]}")
    for name, s in segs.items():
        L = R_OUT * sum(abs(np.angle(np.exp(1j * (np.arctan2(*xyz[b, 1::-1]) - np.arctan2(*xyz[a, 1::-1])))))
                        for a, b in s[:, :2])
        print(f"  {name}: {len(s)} segments, length {L:.5f} m ({np.degrees(L/R_OUT):.2f} deg)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--lc", type=float, default=0.006)
    ap.add_argument("--lc_sheath", type=float, default=0.002)
    ap.add_argument("--lc_center", type=float, default=0.0015)
    ap.add_argument("--order", type=int, default=2, choices=[1, 2])
    ap.add_argument("--w_pee", type=float, default=18.0, help="PEE arc width [deg]")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    build(a.lc, a.lc_sheath, a.lc_center, a.w_pee, a.order, a.out)
