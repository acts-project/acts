#!/usr/bin/env python3
"""Draw an r-z view of a (Gen3) tracking geometry and highlight lossy material merges.

This is a standalone script that only depends on the **pure-core** ACTS Python
bindings (the ``acts`` module, including ``acts.json``) plus matplotlib. It does
*not* use the Examples framework.

It reads a tracking-geometry JSON file (as produced by
``acts.json.TrackingGeometryJsonConverter``), draws every cylinder volume as a
rectangle in the r-z plane, and overlays -- as orange line segments -- any portal
surface that carries a ``MergedMaterialMarker``. The marker is left behind by
``Portal::merge`` when it runs in "keep going" mode and has to discard surface
material during container stacking, so this plot shows exactly *where* those
problematic merges happened.

Example::

    python material_merge_rz_view.py tracking-geometry.json -o detector_rz.svg
"""

import argparse
from pathlib import Path

import acts


class _RZCollector(acts.TrackingGeometryMutableVisitor):
    """Visitor collecting cylinder volume rectangles and marker line segments."""

    def __init__(self, gctx: acts.GeometryContext):
        super().__init__()
        self._gctx = gctx
        # Each entry: (z_min, z_max, r_min, r_max, name)
        self.volumes = []
        # Each entry: (z0, r0, z1, r1)
        self.markers = []
        self._seen_markers = set()

    def visitVolume(self, volume: acts.Volume):
        bounds = volume.volumeBounds
        if bounds.type() != acts.VolumeBoundsType.Cylinder:
            return
        # CylinderVolumeBounds: values()[0..2] = rMin, rMax, halfLengthZ
        values = bounds.values()
        r_min, r_max, half_z = values[0], values[1], values[2]
        z_center = volume.center(self._gctx)[2]
        name = ""
        if isinstance(volume, acts.TrackingVolume):
            name = volume.volumeName
        self.volumes.append((z_center - half_z, z_center + half_z, r_min, r_max, name))

    def visitPortal(self, portal: acts.Portal):
        self._maybe_marker(portal.surface)

    def visitSurface(self, surface: acts.Surface):
        self._maybe_marker(surface)

    def _maybe_marker(self, surface: acts.Surface):
        material = surface.surfaceMaterial
        if not isinstance(material, acts.MergedMaterialMarker):
            return
        segment = self._surface_segment(surface)
        if segment is None:
            return
        # A merged portal surface is reachable both as a portal and (potentially)
        # as a volume surface, and is shared between the stacked volumes, so the
        # same segment is seen multiple times. Deduplicate by rounded geometry.
        key = tuple(round(v, 3) for v in segment)
        if key in self._seen_markers:
            return
        self._seen_markers.add(key)
        self.markers.append(segment)

    def _surface_segment(self, surface: acts.Surface):
        """Return the (z0, r0, z1, r1) r-z footprint of a portal surface."""
        bounds = surface.bounds
        center = surface.center(self._gctx)
        z_center = center[2]
        if isinstance(bounds, acts.CylinderBounds):
            # Cylinder face: horizontal line at fixed radius spanning z.
            r = bounds.values()[int(acts.CylinderBoundsValue.R)]
            half_z = bounds.values()[int(acts.CylinderBoundsValue.HalfLengthZ)]
            return (z_center - half_z, r, z_center + half_z, r)
        if isinstance(bounds, acts.RadialBounds):
            # Disc face: vertical line at fixed z spanning r.
            r_min = bounds.values()[int(acts.RadialBoundsValue.MinR)]
            r_max = bounds.values()[int(acts.RadialBoundsValue.MaxR)]
            return (z_center, r_min, z_center, r_max)
        # Unsupported surface kind for the r-z view; skip it.
        return None


def load_tracking_geometry(
    json_path: Path, gctx: acts.GeometryContext
) -> acts.TrackingGeometry:
    """Load a Gen3 tracking geometry from a JSON file."""
    from acts.json import TrackingGeometryJsonConverter

    converter = TrackingGeometryJsonConverter()
    return converter.fromJson(gctx, json_path.absolute())


def collect_rz(tracking_geometry: acts.TrackingGeometry, gctx: acts.GeometryContext):
    """Collect cylinder volume rectangles and marker segments from the geometry.

    Returns a tuple ``(volumes, markers)``.
    """
    collector = _RZCollector(gctx)
    tracking_geometry.apply(collector)
    return collector.volumes, collector.markers


def render(volumes, markers, output_path: Path, title: str = "Material merge r-z view"):
    """Render the r-z view to an SVG (or any matplotlib-supported) file."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, ax = plt.subplots(figsize=(12, 6))

    for z_min, z_max, r_min, r_max, name in volumes:
        ax.add_patch(
            Rectangle(
                (z_min, r_min),
                z_max - z_min,
                r_max - r_min,
                facecolor="#dfe7f2",
                edgecolor="#34507a",
                linewidth=0.6,
                zorder=1,
            )
        )
        if name:
            # ODD volume names are hierarchical (pipe-separated); show the leaf.
            leaf = name.split("|")[-1]
            # Rotate the label to run along the longer side of the rectangle.
            rotation = 90 if (r_max - r_min) > (z_max - z_min) else 0
            ax.text(
                0.5 * (z_min + z_max),
                0.5 * (r_min + r_max),
                leaf,
                ha="center",
                va="center",
                rotation=rotation,
                fontsize=5,
                color="#1b2a44",
                clip_on=True,
                zorder=3,
            )

    marker_label_used = False
    for z0, r0, z1, r1 in markers:
        ax.plot(
            [z0, z1],
            [r0, r1],
            color="darkorange",
            linewidth=3.0,
            solid_capstyle="round",
            zorder=5,
            label=None if marker_label_used else "Discarded merge material",
        )
        marker_label_used = True

    ax.set_xlabel("z [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title(title)
    ax.autoscale_view()
    ax.margins(0.02)
    ax.set_ylim(bottom=0)
    if marker_label_used:
        ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(str(output_path))
    plt.close(fig)
    return output_path


def run(json_path: Path, output_path: Path, title: str = "Material merge r-z view"):
    """Load a geometry JSON, build the r-z view, and write it to ``output_path``."""
    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    tracking_geometry = load_tracking_geometry(json_path, gctx)
    volumes, markers = collect_rz(tracking_geometry, gctx)
    render(volumes, markers, output_path, title=title)
    return volumes, markers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "geometry",
        type=Path,
        help="Path to a Gen3 tracking-geometry JSON file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("geometry_rz.svg"),
        help="Output image path (default: geometry_rz.svg)",
    )
    parser.add_argument(
        "--title",
        default="Material merge r-z view",
        help="Plot title",
    )
    args = parser.parse_args()

    volumes, markers = run(args.geometry, args.output, title=args.title)
    print(
        f"Drew {len(volumes)} cylinder volumes and {len(markers)} "
        f"merge-material markers to {args.output}"
    )


if __name__ == "__main__":
    main()
