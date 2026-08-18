"""ODD KF alignment workflow — misalignment stage: surface selection, scramble, I/O."""

from __future__ import annotations

import random
from pathlib import Path
from typing import Optional

import acts
import acts.examples.alignment

from kfalignment_config import DEFAULT_MISALIGNMENT_CONFIG, MisalignmentConfig
from kfalignment_geoid import collect_matching_sensitive_surfaces
from kfalignment_transform_io import (
    build_geoid_surface_index_map,
    copy_transform3,
    sort_surfaces_alignment_order,
    write_alignment_delta_res,
    write_misalignment_index_map,
)


def _asymmetric_random(mag: float, enabled: bool) -> float:
    if not enabled:
        return 0.0
    if random.random() < 0.5:
        return random.uniform(-mag * 1.5, -mag * 0.5)
    return random.uniform(mag * 0.5, mag * 1.5)


def apply_random_misalignment_to_transform(
    trf: acts.Transform3,
    tx: int = 1,
    ty: int = 1,
    tz: int = 0,
    rx: int = 0,
    ry: int = 0,
    rz: int = 1,
    shift_mag_mm: float = 0.5,
    rotation_mag_rad: float = 0.02,
) -> tuple:
    """Apply random local tx/ty/rz to ``trf`` in place; return ``(trf, deltas)``."""
    _ = (tz, rx, ry)
    dx = dy = dz = 0.0
    rotx = roty = rotz = 0.0

    if tx != 0:
        l_shift_x = acts.examples.alignment.AlignmentGeneratorLocalShift()
        l_shift_x.axisDirection = acts.AxisDirection.AxisX
        dx = _asymmetric_random(shift_mag_mm, True)
        l_shift_x.shift = dx
        l_shift_x(trf)

    if ty != 0:
        l_shift_y = acts.examples.alignment.AlignmentGeneratorLocalShift()
        l_shift_y.axisDirection = acts.AxisDirection.AxisY
        dy = _asymmetric_random(shift_mag_mm, True)
        l_shift_y.shift = dy
        l_shift_y(trf)

    if rz != 0:
        l_rot_z = acts.examples.alignment.AlignmentGeneratorLocalRotation()
        l_rot_z.axis = acts.Vector3(0.0, 0.0, 1.0)
        rotz = _asymmetric_random(rotation_mag_rad, True)
        l_rot_z.angle = rotz
        l_rot_z(trf)

    return trf, [dx, dy, dz, rotx, roty, rotz]


def _print_misalignment_filter_criteria(
    target_volume,
    target_layer,
    target_sensitive,
    target_extra,
    num_selected: int,
) -> None:
    print(
        f"Misalignment filter: volume={target_volume!r}, "
        f"layer={target_layer}, sensitive={target_sensitive}, extra={target_extra} "
        f"-> {num_selected} surfaces"
    )


def setupMisalignment(
    trackingGeometry: acts.TrackingGeometry,
    outputDir: Path,
    config: Optional[MisalignmentConfig] = None,
    target_volume=-1,
    target_layer=-1,
    target_sensitive=-1,
    target_extra: int = -1,
    tx: int = 1,
    ty: int = 1,
    tz: int = 0,
    rx: int = 0,
    ry: int = 0,
    rz: int = 1,
    shift_mag_mm: float = 0.5,
    rotation_mag_rad: float = 0.02,
    seed: int = 42,
):
    """Apply random misalignment and write applied + index text files."""
    if config is not None:
        target_volume = config.target_volume
        target_layer = config.target_layer
        target_sensitive = config.target_sensitive
        target_extra = config.target_extra
        tx, ty, tz = config.tx, config.ty, config.tz
        rx, ry, rz = config.rx, config.ry, config.rz
        shift_mag_mm = config.shift_mag_mm
        rotation_mag_rad = config.rotation_mag_rad
        seed = config.seed

    random.seed(seed)

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    matches = collect_matching_sensitive_surfaces(
        trackingGeometry,
        target_volume,
        target_layer,
        target_sensitive,
        target_extra,
    )
    matches = sort_surfaces_alignment_order(matches)
    index_map = build_geoid_surface_index_map(matches)

    _print_misalignment_filter_criteria(
        target_volume,
        target_layer,
        target_sensitive,
        target_extra,
        len(matches),
    )
    print(f"Configured DoF mask: tx={tx} ty={ty} tz={tz} " f"rx={rx} ry={ry} rz={rz}")
    print(
        f"Random scramble DoF (only tx/ty/rz): tx={int(bool(tx))} "
        f"ty={int(bool(ty))} rz={int(bool(rz))}, "
        f"shift=±{shift_mag_mm} mm, rotation=±{rotation_mag_rad} rad, "
        f"seed={seed}"
    )

    geoIdMap = {}
    alignment_placements = []
    surface_deltas = []

    for surface, placement in matches:
        gid = surface.geometryId
        trf = copy_transform3(surface.localToGlobalTransform(gctx))
        trf, deltas = apply_random_misalignment_to_transform(
            trf,
            tx=tx,
            ty=ty,
            tz=0,
            rx=0,
            ry=0,
            rz=rz,
            shift_mag_mm=shift_mag_mm,
            rotation_mag_rad=rotation_mag_rad,
        )
        geoIdMap[gid] = trf
        alignment_placements.append(placement)
        surface_deltas.append(deltas)

    outputDir = Path(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)

    applied_file = outputDir / "misalignment_applied.txt"
    index_file = outputDir / "misalignment_index_map.txt"
    write_alignment_delta_res(applied_file, surface_deltas)
    write_misalignment_index_map(index_file, index_map)

    print(
        f"Misaligned {len(geoIdMap)} elements; wrote "
        f"{applied_file.name}, {index_file.name}"
    )

    return geoIdMap, alignment_placements


def runMisalignment(
    trackingGeometry: acts.TrackingGeometry,
    outputDir: Path,
    config: Optional[MisalignmentConfig] = None,
):
    """CLI/full-workflow entry: write misalignment files under ``outputDir``."""
    cfg = config if config is not None else DEFAULT_MISALIGNMENT_CONFIG
    outputDir = Path(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)
    print(f"Misalignment-only: writing outputs to {outputDir}")
    return setupMisalignment(
        trackingGeometry=trackingGeometry,
        outputDir=outputDir,
        config=cfg,
    )
