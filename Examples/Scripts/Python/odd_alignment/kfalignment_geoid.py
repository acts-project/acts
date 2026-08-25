"""Geometry-ID surface selection for misalignment / alignment."""

from __future__ import annotations

import acts
import acts.examples.alignment


def match_surface_to_filter(gid, target_vol, target_lay, target_sen, target_ext):
    """Return True if ``gid`` matches the misalignment surface filter."""
    vol_id = gid.volume
    lay_id = gid.layer
    sen_id = gid.sensitive

    if isinstance(target_vol, dict):
        if vol_id not in target_vol:
            return False

        layer_config = target_vol[vol_id]

        if layer_config == -1:
            pass
        elif isinstance(layer_config, dict):
            if lay_id not in layer_config:
                return False

            sensor_config = layer_config[lay_id]

            if sensor_config == -1:
                pass
            elif isinstance(sensor_config, list):
                if sen_id not in sensor_config:
                    return False
            else:
                if sen_id != sensor_config:
                    return False
        elif isinstance(layer_config, list):
            if lay_id not in layer_config:
                return False
        else:
            if lay_id != layer_config:
                return False
    else:
        if target_vol != -1:
            if isinstance(target_vol, list):
                if vol_id not in target_vol:
                    return False
            else:
                if vol_id != target_vol:
                    return False

        if target_lay != -1:
            if isinstance(target_lay, list):
                if lay_id not in target_lay:
                    return False
            else:
                if lay_id != target_lay:
                    return False

        if target_sen != -1:
            if isinstance(target_sen, list):
                if sen_id not in target_sen:
                    return False
            else:
                if sen_id != target_sen:
                    return False

    if target_ext != -1:
        if gid.extra != target_ext:
            return False
    return True


def collect_matching_sensitive_surfaces(
    trackingGeometry: acts.TrackingGeometry,
    target_volume=-1,
    target_layer=-1,
    target_sensitive=-1,
    target_extra: int = -1,
) -> list:
    """Collect (surface, placement) pairs matching the filter (single geometry visit)."""
    matches = []

    def visit(surface: acts.Surface) -> bool:
        if not surface.isSensitive:
            return True
        placement = acts.examples.alignment.surfacePlacement(surface)
        if placement is None:
            return True
        if match_surface_to_filter(
            surface.geometryId,
            target_volume,
            target_layer,
            target_sensitive,
            target_extra,
        ):
            matches.append((surface, placement))
        return True

    trackingGeometry.visitSurfaces(visit)
    return matches
