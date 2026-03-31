"""2D projections of surface polyhedra (xy, zr, zphi)."""

from __future__ import annotations

import math
from typing import Iterable, List, Literal, Sequence, Tuple

import numpy as np

ViewName = Literal["xy", "zr", "zphi"]


def project_point(
    view: ViewName,
    x: float,
    y: float,
    z: float,
) -> Tuple[float, float]:
    if view == "xy":
        return x, y
    if view == "zr":
        r = math.hypot(x, y)
        return z, r
    if view == "zphi":
        phi = math.atan2(y, x)
        return z, phi
    raise ValueError(view)


def project_polyhedron_vertices(
    view: ViewName,
    vertices,
) -> np.ndarray:
    """Project ACTS ``Polyhedron.vertices`` (Vector3 sequence) to Nx2 float array."""

    pts = np.zeros((len(vertices), 2))
    for i, v in enumerate(vertices):
        pts[i] = project_point(view, float(v[0]), float(v[1]), float(v[2]))
    return pts


def _dedupe_loop(idxs: Sequence[int]) -> List[int]:
    if not idxs:
        return []
    out = list(idxs)
    if out[0] == out[-1]:
        out = out[:-1]
    return out


def surface_to_2d_patches(
    view: ViewName,
    polyhedron,
    *,
    prefer_triangles: bool = False,
) -> List[np.ndarray]:
    """Return closed 2D polygons (each Kx2) for matplotlib.

    By default uses the native ``Polyhedron.faces`` (polygon mesh). With
    ``prefer_triangles=True``, uses ``triangularMesh`` when non-empty.
    """

    if prefer_triangles and polyhedron.triangularMesh:
        faces_source = polyhedron.triangularMesh
    else:
        faces_source = polyhedron.faces
    patches: List[np.ndarray] = []
    verts3 = polyhedron.vertices

    for face in faces_source:
        idxs = _dedupe_loop(face)
        ring = np.zeros((len(idxs), 2))
        for j, k in enumerate(idxs):
            v = verts3[k]
            ring[j] = project_point(view, float(v[0]), float(v[1]), float(v[2]))
        if len(ring) >= 3:
            patches.append(ring)
    return patches


def auto_axis_limits(
    polygons: Iterable[np.ndarray],
    *,
    margin_ratio: float = 0.05,
) -> Tuple[float, float, float, float]:
    """Return xmin, ymin, xmax, ymax over all polygon vertices."""

    xs: List[float] = []
    ys: List[float] = []
    for poly in polygons:
        if poly.size == 0:
            continue
        xs.extend(poly[:, 0].tolist())
        ys.extend(poly[:, 1].tolist())
    if not xs:
        return -1.0, -1.0, 1.0, 1.0
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    dx = xmax - xmin if xmax > xmin else 1.0
    dy = ymax - ymin if ymax > ymin else 1.0
    mx = margin_ratio * dx
    my = margin_ratio * dy
    return xmin - mx, ymin - my, xmax + mx, ymax + my
