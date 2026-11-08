"""
Circular 2D triangle mesh for an OGS T (heat transport) model - STIMTEC.

Geometry
  - outer circle  R_out = 0.100 m  (model boundary)
  - inner circle  R_in  = 0.066 m  (embedded; used later for volumetric strain)
  - central node  (0, 0, 0)        (embedded point, always a mesh node)

Output (VTU, OGS-ready)
  domain.vtu         bulk mesh, linear triangles, MaterialIDs = 0
  outer_boundary.vtu line elements on r = R_out
  inner_circle.vtu   line elements on r = R_in
  inner_disk.vtu     triangles with r <= R_in (for averaging/integration)
  center.vtu         single node at (0,0,0)
All subdomain meshes carry bulk_node_ids (and bulk_element_ids for
line/area meshes), so they can be used directly as OGS <meshes>.

Usage:  python make_circle_mesh.py [--lc 0.004] [--lc_center 0.002] [--out .]
"""
import argparse
import os

import gmsh
import meshio
import numpy as np

R_OUT = 0.1
R_IN = 0.066


def build(lc, lc_center, outdir):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("stimtec_circle")
    occ = gmsh.model.occ

    outer = occ.addDisk(0, 0, 0, R_OUT, R_OUT)
    inner = occ.addDisk(0, 0, 0, R_IN, R_IN)
    # fragment -> inner disk + annulus sharing the r = R_in circle (conforming)
    occ.fragment([(2, outer)], [(2, inner)])
    pc = occ.addPoint(0, 0, 0, lc_center)
    occ.synchronize()

    # identify surfaces / curves by geometry
    surfs = gmsh.model.getEntities(2)
    def surf_rmax(tag):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(2, tag)
        return max(xmax, ymax)
    s_inner = [t for _, t in surfs if abs(surf_rmax(t) - R_IN) < 1e-4]
    s_ann = [t for _, t in surfs if abs(surf_rmax(t) - R_OUT) < 1e-4]
    curves = gmsh.model.getEntities(1)
    def curve_r(tag):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(1, tag)
        return max(xmax, ymax)
    c_out = [t for _, t in curves if abs(curve_r(t) - R_OUT) < 1e-4]
    c_in = [t for _, t in curves if abs(curve_r(t) - R_IN) < 1e-4]

    # embed centre point so (0,0,0) is a mesh node
    gmsh.model.mesh.embed(0, [pc], 2, s_inner[0])

    # size field: lc_center at centre, growing linearly to lc at R_IN
    f = gmsh.model.mesh.field
    fd = f.add("Distance"); f.setNumbers(fd, "PointsList", [pc])
    ft = f.add("Threshold"); f.setNumber(ft, "InField", fd)
    f.setNumber(ft, "SizeMin", lc_center); f.setNumber(ft, "SizeMax", lc)
    f.setNumber(ft, "DistMin", 0.0); f.setNumber(ft, "DistMax", R_IN)
    f.setAsBackgroundMesh(ft)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay

    gmsh.model.addPhysicalGroup(2, s_inner, 1, "inner_disk")
    gmsh.model.addPhysicalGroup(2, s_ann, 2, "annulus")
    gmsh.model.addPhysicalGroup(1, c_out, 3, "outer_boundary")
    gmsh.model.addPhysicalGroup(1, c_in, 4, "inner_circle")
    gmsh.model.addPhysicalGroup(0, [pc], 5, "center")

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Laplace2D")
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()
    gmsh.write(os.path.join(outdir, "circle.msh"))

    # ---------------- nodes ----------------
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    order = np.argsort(tags)
    tags, xyz = tags[order], xyz[order]
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

    tri_in = elems(2, 1)
    tri_an = elems(2, 2)
    seg_out = elems(1, 3)
    seg_in = elems(1, 4)
    center_node = int(t2i[gmsh.model.mesh.getNodesForPhysicalGroup(0, 5)[0][0]])
    gmsh.finalize()

    # ensure counter-clockwise triangles (positive area, required by OGS)
    tri = np.vstack([tri_in, tri_an])
    p = xyz[tri]
    area2 = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - \
            (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])
    tri[area2 < 0] = tri[area2 < 0][:, [0, 2, 1]]
    n_in = len(tri_in)

    # ---------------- bulk mesh ----------------
    bulk = meshio.Mesh(
        xyz, [("triangle", tri)],
        cell_data={"MaterialIDs": [np.zeros(len(tri), dtype=np.int32)]})
    meshio.write(os.path.join(outdir, "domain.vtu"), bulk, binary=True)

    # ---------------- subdomain meshes ----------------
    # map boundary segment -> adjacent bulk triangle (for bulk_element_ids)
    edge2tri = {}
    for e, t in enumerate(tri):
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            edge2tri.setdefault(frozenset((a, b)), e)

    def write_sub(name, cells, celltype, bulk_eids=None):
        used = np.unique(cells)
        loc = {g: i for i, g in enumerate(used)}
        lc_ = np.vectorize(loc.get)(cells)
        pd = {"bulk_node_ids": used.astype(np.uint64)}
        cd = {}
        if bulk_eids is not None:
            cd["bulk_element_ids"] = [np.asarray(bulk_eids, dtype=np.uint64)]
        m = meshio.Mesh(xyz[used], [(celltype, lc_)], point_data=pd, cell_data=cd)
        meshio.write(os.path.join(outdir, name), m, binary=True)

    write_sub("outer_boundary.vtu", seg_out, "line",
              [edge2tri[frozenset(s)] for s in seg_out])
    write_sub("inner_circle.vtu", seg_in, "line",
              [edge2tri[frozenset(s)] for s in seg_in])
    write_sub("inner_disk.vtu", tri[:n_in], "triangle", np.arange(n_in))
    write_sub("center.vtu", np.array([[center_node]]), "vertex")

    # ---------------- report ----------------
    r = np.hypot(xyz[:, 0], xyz[:, 1])
    A = 0.5 * np.abs(area2)
    print(f"nodes {len(xyz)}, triangles {len(tri)} (inner {n_in}, annulus {len(tri)-n_in})")
    print(f"outer segments {len(seg_out)}, inner-circle segments {len(seg_in)}")
    print(f"centre node id {center_node}: {xyz[center_node]}")
    print(f"max r = {r.max():.6f}, nodes on R_in: {np.sum(np.abs(r-R_IN)<1e-9)}")
    print(f"area mesh {A.sum():.6e}  vs pi R^2 {np.pi*R_OUT**2:.6e}")
    print(f"area inner {A[:n_in].sum():.6e} vs pi r^2 {np.pi*R_IN**2:.6e}")
    return xyz, tri, n_in, seg_out, seg_in, center_node


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--lc", type=float, default=0.004, help="element size at r >= R_in [m]")
    ap.add_argument("--lc_center", type=float, default=0.002, help="element size at centre [m]")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    build(a.lc, a.lc_center, a.out)
