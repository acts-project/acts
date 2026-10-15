#!/usr/bin/env python3
"""Export actual gen1 ODD material surfaces for the core-only C++ benchmark.

Run in an ACTS Python environment with DD4hep, ROOT and ODD_PATH configured.
The deliberately small text format stores geometry IDs, all four bounds values,
and the complete placement (row-major 3x3 rotation followed by translation).
"""

import argparse
from collections import Counter

import acts
from acts.examples.odd import getOpenDataDetector

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output")
args = parser.parse_args()
detector = getOpenDataDetector(gen3=False, logLevel=acts.logging.ERROR)
surfaces = detector.trackingGeometry().extractMaterialSurfaces()
gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
with open(args.output, "w") as output:
    output.write(
        "# gen1 ODD extractMaterialSurfaces; mm; kind id bounds[4] rotation[9] translation[3]\n"
    )
    for surface in surfaces:
        if surface.type == acts.SurfaceType.Cylinder:
            kind = "cylinder"
        elif (
            surface.type == acts.SurfaceType.Disc
            and surface.bounds.type == acts.SurfaceBoundsType.Disc
        ):
            kind = "disc"
        else:
            raise RuntimeError(f"Unsupported surface: {surface.type}, {surface.bounds}")
        transform = surface.localToGlobalTransform(gctx)
        values = [
            *surface.bounds.values(),
            *(transform.rotation[i, j] for i in range(3) for j in range(3)),
            *transform.translation,
        ]
        output.write(
            f"{kind} {surface.geometryId.value} "
            + " ".join(format(float(v), ".17g") for v in values)
            + "\n"
        )
print(f"Exported {len(surfaces)} surfaces: {Counter(str(s.type) for s in surfaces)}")
