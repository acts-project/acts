from pathlib import Path
from typing import Optional
import json
import random
import numpy as np
import ctypes
import os
import argparse
import sys


def _preload_dd4hep_for_macos() -> None:
    """Load DD4hep Gaudi plugin + core plugins with RTLD_GLOBAL before acts.

    On macOS, Python-bundled .so and DD4hep's Gaudi plugin registry can get out
    of sync (empty stub / bad any_cast for lccdd_XML_reader). Preloading the
    plugin manager and libDDCorePlugins in the same process often fixes that.
    """
    if sys.platform != "darwin":
        return
    root = os.environ.get("DD4HEP_ROOT") or os.environ.get("DD4HEP_INSTALL")
    libdir: Optional[Path] = Path(root) / "lib" if root else None
    if libdir is None or not libdir.is_dir():
        try:
            candidate = (
                Path(__file__).resolve().parents[5] / "DD4hep" / "install" / "lib"
            )
        except (IndexError, OSError):
            candidate = None
        if candidate is not None and candidate.is_dir():
            libdir = candidate
    if libdir is None or not libdir.is_dir():
        return
    path_bits = [str(libdir)]
    odd = os.environ.get("ODD_PATH")
    if odd:
        fdir = Path(odd).resolve().parent / "odd-build" / "factory"
        if fdir.is_dir():
            path_bits.append(str(fdir))
    for var in ("DYLD_LIBRARY_PATH", "DD4HEP_LIBRARY_PATH"):
        cur = os.environ.get(var, "")
        merged = []
        for p in path_bits + [x for x in cur.split(":") if x]:
            if p not in merged:
                merged.append(p)
        os.environ[var] = ":".join(merged)
    mode = getattr(ctypes, "RTLD_GLOBAL", 8)
    for names in (
        ("libDD4hepGaudiPluginMgr.1.32.dylib", "libDD4hepGaudiPluginMgr.dylib"),
        ("libDDCorePlugins.1.32.dylib", "libDDCorePlugins.dylib"),
    ):
        for n in names:
            p = libdir / n
            if p.is_file():
                try:
                    ctypes.CDLL(str(p), mode=mode)
                except OSError:
                    break
                break


_preload_dd4hep_for_macos()

import acts
import acts.examples
import acts.examples.alignment
from acts.examples.alignment import AlignmentDecorator

# Sequencer may otherwise raise after Kalman+RootTrackStatesWriter (FLTINV in writeT).
# ACTS CI uses ACTS_SEQUENCER_FAIL_ON_UNMASKED_FPE=0; set to 1 to fail on any unmasked FPE.
if os.environ.get("ACTS_SEQUENCER_FAIL_ON_UNMASKED_FPE") is None:
    os.environ["ACTS_SEQUENCER_FAIL_ON_UNMASKED_FPE"] = "0"

u = acts.UnitConstants


def parse_geoid(id_str):
    """Parse ID string, e.g. 'vol=16|lay=4|sen=1|ext=2'"""
    parts = {}
    for part in id_str.split("|"):
        key, value = part.split("=")
        parts[key] = int(value)
    return parts


def geoid_to_id_string(gid: acts.GeometryIdentifier) -> str:
    """Convert GeometryIdentifier to ID string format: 'vol=X|lay=Y|sen=Z|ext=W'"""
    parts = [f"vol={gid.volume}", f"lay={gid.layer}", f"sen={gid.sensitive}"]
    if gid.extra != 0:
        parts.append(f"ext={gid.extra}")
    return "|".join(parts)


def match_filter(geoid_parts, target_vol, target_lay, target_sen, target_ext):
    """
    Check if geoid matches filtering conditions

    Note: Some elements don't have 'ext' parameter (e.g., volume=17,24,29)
    - If target_ext == -1: match all elements (with or without ext)
    - If target_ext != -1: only match elements that have ext AND ext == target_ext

    target_vol can be:
    - -1: match all volumes
    - int: match specific volume
    - list: match any volume in the list
    - dict: hierarchical structure {volume: {layer: [sensors]}}
        Example: {16: {4: [1, 2, 3], 5: [1, 2]}, 17: {2: [1]}}
        If layer value is -1 or list contains -1, match all sensors in that layer
        If layer value is a list, match sensors in any of those layers

    target_lay, target_sen: used for backward compatibility when target_vol is not a dict
    - -1: match all
    - int: match specific value
    - list: match any value in the list
    """
    vol_id = geoid_parts.get("vol")
    lay_id = geoid_parts.get("lay")
    sen_id = geoid_parts.get("sen")

    # Support hierarchical structure: {volume: {layer: [sensors]}}
    if isinstance(target_vol, dict):
        # Check if volume exists in config
        if vol_id not in target_vol:
            return False

        layer_config = target_vol[vol_id]

        # If layer_config is -1, match all layers and sensors in this volume
        if layer_config == -1:
            pass  # Continue to check ext
        elif isinstance(layer_config, dict):
            # Check if layer exists in config
            if lay_id not in layer_config:
                return False

            sensor_config = layer_config[lay_id]

            # If sensor_config is -1, match all sensors in this layer
            if sensor_config == -1:
                pass  # Continue to check ext
            elif isinstance(sensor_config, list):
                # Check if sensor is in the list
                if sen_id not in sensor_config:
                    return False
            else:
                # Single sensor value
                if sen_id != sensor_config:
                    return False
        elif isinstance(layer_config, list):
            # List of layers: match sensors in any of these layers
            if lay_id not in layer_config:
                return False
        else:
            # Single layer value
            if lay_id != layer_config:
                return False
    else:
        # Backward compatibility: original flat structure
        # Support target_vol as list or single value
        if target_vol != -1:
            if isinstance(target_vol, list):
                if vol_id not in target_vol:
                    return False
            else:
                if vol_id != target_vol:
                    return False

        # Support target_lay as list or single value
        if target_lay != -1:
            if isinstance(target_lay, list):
                if lay_id not in target_lay:
                    return False
            else:
                if lay_id != target_lay:
                    return False

        # Support target_sen as list or single value
        if target_sen != -1:
            if isinstance(target_sen, list):
                if sen_id not in target_sen:
                    return False
            else:
                if sen_id != target_sen:
                    return False

    if target_ext != -1:
        if "ext" not in geoid_parts or geoid_parts["ext"] != target_ext:
            return False
    return True


def match_surface_to_filter(gid, target_vol, target_lay, target_sen, target_ext):
    """
    Check if surface's GeoID matches filtering conditions

    target_vol can be:
    - -1: match all volumes
    - int: match specific volume
    - list: match any volume in the list
    - dict: hierarchical structure {volume: {layer: [sensors]}}
        Example: {16: {4: [1, 2, 3], 5: [1, 2]}, 17: {2: [1]}}
        If layer value is -1 or list contains -1, match all sensors in that layer
        If layer value is a list, match sensors in any of those layers

    target_lay, target_sen: used for backward compatibility when target_vol is not a dict
    - -1: match all
    - int: match specific value
    - list: match any value in the list
    """
    vol_id = gid.volume
    lay_id = gid.layer
    sen_id = gid.sensitive

    # Support hierarchical structure: {volume: {layer: [sensors]}}
    if isinstance(target_vol, dict):
        # Check if volume exists in config
        if vol_id not in target_vol:
            return False

        layer_config = target_vol[vol_id]

        # If layer_config is -1, match all layers and sensors in this volume
        if layer_config == -1:
            pass  # Continue to check ext
        elif isinstance(layer_config, dict):
            # Check if layer exists in config
            if lay_id not in layer_config:
                return False

            sensor_config = layer_config[lay_id]

            # If sensor_config is -1, match all sensors in this layer
            if sensor_config == -1:
                pass  # Continue to check ext
            elif isinstance(sensor_config, list):
                # Check if sensor is in the list
                if sen_id not in sensor_config:
                    return False
            else:
                # Single sensor value
                if sen_id != sensor_config:
                    return False
        elif isinstance(layer_config, list):
            # List of layers: match sensors in any of these layers
            if lay_id not in layer_config:
                return False
        else:
            # Single layer value
            if lay_id != layer_config:
                return False
    else:
        # Backward compatibility: original flat structure
        # Support target_vol as list or single value
        if target_vol != -1:
            if isinstance(target_vol, list):
                if vol_id not in target_vol:
                    return False
            else:
                if vol_id != target_vol:
                    return False

        # Support target_lay as list or single value
        if target_lay != -1:
            if isinstance(target_lay, list):
                if lay_id not in target_lay:
                    return False
            else:
                if lay_id != target_lay:
                    return False

        # Support target_sen as list or single value
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


# ============================================================================
# Simulation functions
# ============================================================================


def runSimulation(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    detector=None,
    numEvents=100,
    s: acts.examples.Sequencer = None,
):
    """
    Run simulation with misaligned geometry.
    """
    from acts.examples.simulation import (
        addParticleGun,
        ParticleConfig,
        EtaConfig,
        PhiConfig,
        MomentumConfig,
        addFatras,
        addDigitization,
    )

    s = s or acts.examples.Sequencer(
        events=numEvents, numThreads=-1, logLevel=acts.logging.INFO
    )

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)
    outputDir.mkdir(exist_ok=True)

    logger = acts.getDefaultLogger("Simulation", acts.logging.INFO)

    addParticleGun(
        s,
        ParticleConfig(num=1, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
        EtaConfig(-3.0, 3.0, uniform=True),
        MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
        PhiConfig(0.0 * u.degree, 360.0 * u.degree),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(2.0 * u.mm, 2.0 * u.mm, 2.0 * u.mm, 2.0 * u.ns),
        ),
        multiplicity=1,
        rnd=rnd,
        outputDirCsv=outputDir,
        outputDirRoot=outputDir,
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        enableInteractions=True,
        outputDirCsv=outputDir,
        outputDirRoot=outputDir,
    )

    digiAlg = addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
        outputDirCsv=outputDir,
        outputDirRoot=outputDir,
    )

    # Note: particle_measurements_map will be generated after simulation
    # using the post-processing function below
    logger.info("Simulation output will be saved to: {}", str(outputDir))
    logger.info(
        "particle_measurements_map will be generated after simulation completes"
    )

    return s


# ============================================================================
# Alignment functions
# ============================================================================


def copy_transform3(trf: acts.Transform3) -> acts.Transform3:
    """Return a mutable copy of an Acts Transform3."""
    t = trf.translation
    r = trf.rotation
    return acts.Transform3(
        acts.Vector3(float(t[0]), float(t[1]), float(t[2])),
        acts.RotationMatrix3(
            acts.Vector3(r[0, 0], r[1, 0], r[2, 0]),
            acts.Vector3(r[0, 1], r[1, 1], r[2, 1]),
            acts.Vector3(r[0, 2], r[1, 2], r[2, 2]),
        ),
    )


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


def _asymmetric_random(mag: float, enabled: bool) -> float:
    if not enabled:
        return 0.0
    if random.random() < 0.5:
        return random.uniform(-mag * 1.5, -mag * 0.5)
    return random.uniform(mag * 0.5, mag * 1.5)


def apply_random_misalignment_to_transform(
    trf: acts.Transform3,
    tx: int = 0,
    ty: int = 0,
    tz: int = 0,
    rx: int = 0,
    ry: int = 0,
    rz: int = 0,
    shift_mag_mm: float = 2.0,
    rotation_mag_rad: float = 0.02,
) -> acts.Transform3:
    """Apply random local shifts/rotations to a transform (modified in place)."""
    if tx != 0:
        l_shift_x = acts.examples.alignment.AlignmentGeneratorLocalShift()
        l_shift_x.axisDirection = acts.AxisDirection.AxisX
        l_shift_x.shift = _asymmetric_random(shift_mag_mm, True)
        l_shift_x(trf)

    if ty != 0:
        l_shift_y = acts.examples.alignment.AlignmentGeneratorLocalShift()
        l_shift_y.axisDirection = acts.AxisDirection.AxisY
        l_shift_y.shift = _asymmetric_random(shift_mag_mm, True)
        l_shift_y(trf)

    if tz != 0:
        l_shift_z = acts.examples.alignment.AlignmentGeneratorLocalShift()
        l_shift_z.axisDirection = acts.AxisDirection.AxisZ
        l_shift_z.shift = _asymmetric_random(shift_mag_mm, True)
        l_shift_z(trf)

    current_rotation = trf.rotation

    if rx != 0:
        l_rot_x = acts.examples.alignment.AlignmentGeneratorLocalRotation()
        local_x_axis = acts.Vector3(
            current_rotation[0, 0], current_rotation[1, 0], current_rotation[2, 0]
        )
        l_rot_x.axis = local_x_axis
        l_rot_x.angle = _asymmetric_random(rotation_mag_rad, True)
        l_rot_x(trf)
        current_rotation = trf.rotation

    if ry != 0:
        l_rot_y = acts.examples.alignment.AlignmentGeneratorLocalRotation()
        local_y_axis = acts.Vector3(
            current_rotation[0, 1], current_rotation[1, 1], current_rotation[2, 1]
        )
        l_rot_y.axis = local_y_axis
        l_rot_y.angle = _asymmetric_random(rotation_mag_rad, True)
        l_rot_y(trf)
        current_rotation = trf.rotation

    if rz != 0:
        l_rot_z = acts.examples.alignment.AlignmentGeneratorLocalRotation()
        local_z_axis = acts.Vector3(
            current_rotation[0, 2], current_rotation[1, 2], current_rotation[2, 2]
        )
        l_rot_z.axis = local_z_axis
        l_rot_z.angle = _asymmetric_random(rotation_mag_rad, True)
        l_rot_z(trf)

    return trf


def transform_to_alignment_record(
    gid: acts.GeometryIdentifier, trf: acts.Transform3
) -> dict:
    """Build a JSON-serializable transform record for misalignment/alignment I/O."""
    rot = trf.rotation
    trans = trf.translation
    rotation_matrix = [[float(rot[i, j]) for j in range(3)] for i in range(3)]
    return {
        "ID": geoid_to_id_string(gid),
        "Translation": [float(trans[0]), float(trans[1]), float(trans[2])],
        "RotationMatrix": rotation_matrix,
    }


def geo_id_map_from_transform_records(records: list) -> dict:
    """Load GeometryIdentifier -> Transform3 map from alignment/misalignment records."""
    geo_id_map = {}
    for element in records:
        geoid_parts = parse_geoid(element["ID"])
        gid_kwargs = {
            "volume": geoid_parts.get("vol", 0),
            "layer": geoid_parts.get("lay", 0),
            "sensitive": geoid_parts.get("sen", 0),
        }
        if "ext" in geoid_parts:
            gid_kwargs["extra"] = geoid_parts["ext"]
        gid = acts.GeometryIdentifier(**gid_kwargs)

        trans_data = element["Translation"]
        rot_matrix = element["RotationMatrix"]
        trans = acts.Vector3(trans_data[0], trans_data[1], trans_data[2])
        R = np.array(rot_matrix)
        rotation = acts.RotationMatrix3(
            acts.Vector3(R[0, 0], R[1, 0], R[2, 0]),
            acts.Vector3(R[0, 1], R[1, 1], R[2, 1]),
            acts.Vector3(R[0, 2], R[1, 2], R[2, 2]),
        )
        geo_id_map[gid] = acts.Transform3(trans, rotation)
    return geo_id_map


def _print_misalignment_filter_criteria(
    target_volume,
    target_layer,
    target_sensitive,
    target_extra,
    num_selected: int,
) -> None:
    """Print a one-line misalignment filter summary."""
    print(
        f"Misalignment filter: volume={target_volume!r}, "
        f"layer={target_layer}, sensitive={target_sensitive}, extra={target_extra} "
        f"-> {num_selected} surfaces"
    )


def setupMisalignment(
    trackingGeometry: acts.TrackingGeometry,
    outputDir: Path,
    target_volume=-1,  # Can be int or list of ints
    target_layer=-1,  # Can be int or list of ints
    target_sensitive=-1,  # Can be int or list of ints
    target_extra: int = -1,
    tx: int = 0,
    ty: int = 0,
    tz: int = 1,
    rx: int = 0,
    ry: int = 0,
    rz: int = 0,
    shift_mag_mm: float = 2.0,
    rotation_mag_rad: float = 0.02,
):
    """
    Setup misalignment for selected sensitive surfaces.

    Nominal transforms are read from ``surface.localToGlobalTransform`` on the
    constructed ODD tracking geometry (no odd_transforms.json input).

    Returns
    -------
    tuple
        (geoIdMap, alignment_placements, misalignment_record)
    """
    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    matches = collect_matching_sensitive_surfaces(
        trackingGeometry,
        target_volume,
        target_layer,
        target_sensitive,
        target_extra,
    )
    _print_misalignment_filter_criteria(
        target_volume,
        target_layer,
        target_sensitive,
        target_extra,
        len(matches),
    )
    print(
        f"Misalignment DOF: tx={tx} ty={ty} tz={tz} rx={rx} ry={ry} rz={rz}, "
        f"shift=±{shift_mag_mm} mm, rotation=±{rotation_mag_rad} rad"
    )

    geoIdMap = {}
    misalignment_record = []
    alignment_placements = []

    for surface, placement in matches:
        gid = surface.geometryId
        trf = copy_transform3(surface.localToGlobalTransform(gctx))
        apply_random_misalignment_to_transform(
            trf,
            tx=tx,
            ty=ty,
            tz=tz,
            rx=rx,
            ry=ry,
            rz=rz,
            shift_mag_mm=shift_mag_mm,
            rotation_mag_rad=rotation_mag_rad,
        )
        geoIdMap[gid] = trf
        alignment_placements.append(placement)
        misalignment_record.append(transform_to_alignment_record(gid, trf))

    misalignment_file = outputDir / "misalignment_applied.json"
    with open(misalignment_file, "w") as f:
        json.dump(misalignment_record, f, indent=2)
    print(f"Misaligned {len(geoIdMap)} elements; wrote {misalignment_file.name}")

    return geoIdMap, alignment_placements, misalignment_record


def runAlignment(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    inputDir: Path,
    outputDir: Path,
    detector=None,
    numEvents=100,
    reverseFilteringMomThreshold=0 * u.GeV,
    reverseFilteringCovarianceScaling=1,
    s: acts.examples.Sequencer = None,
):
    from acts.examples.simulation import (
        addFatras,
        addDigitization,
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        addKalmanTracks,
    )

    s = s or acts.examples.Sequencer(
        events=numEvents, numThreads=1, logLevel=acts.logging.INFO
    )

    rnd = acts.examples.RandomNumbers(seed=42)
    inputDir = Path(inputDir)
    outputDir = Path(outputDir)

    logger = acts.getDefaultLogger("Alignment", acts.logging.INFO)

    # Setup misalignment parameters (nominal transforms from tracking geometry)
    # Format 1: Hierarchical structure - specify volume -> layer -> sensors
    # Example: {16: {4: [1, 2, 3], 5: [1, 2]}, 17: {2: [1]}}
    #   - Volume 16, Layer 4: sensors [1, 2, 3]
    #   - Volume 16, Layer 5: sensors [1, 2]
    #   - Volume 17, Layer 2: sensor [1]
    # Use -1 for layer value to match all sensors in that layer
    # Use -1 for sensor list to match all sensors: {16: {4: -1}} means all sensors in vol16/lay4
    # Format 1: Use hierarchical structure (dict) - specify volume -> layer -> sensors
    target_volume = {17: {8: -1}, 24: {4: -1}}
    # Format 2: Use flat structure (backward compatible)
    # target_volume = -1      # -1 = all volumes, or dict for hierarchical: {vol: {lay: [sensors]}}
    # target_layer = -1        # -1 = all layers (used only if target_volume is not a dict)
    # target_sensitive = -1   # Used only if target_volume is not a dict
    # target_extra = -1       # -1 = all extra

    # These variables must be defined even when using hierarchical structure (for backward compatibility)
    target_layer = -1  # -1 = all layers (used only if target_volume is not a dict)
    target_sensitive = -1  # Used only if target_volume is not a dict
    target_extra = -1  # -1 = all extra

    initialvarinflation = [100, 100, 100, 100, 100, 1]

    tx, ty, tz = 1, 1, 0  # Translation DOF
    rx, ry, rz = 1, 1, 1  # Rotation DOF
    shift_mag_mm = 0.5  # Shift magnitude (mm)
    rotation_mag_rad = 0.02  # Rotation magnitude (rad)

    geoIdMap, alignment_placements, _misalignment_record = setupMisalignment(
        trackingGeometry=trackingGeometry,
        outputDir=outputDir,
        target_volume=target_volume,
        target_layer=target_layer,
        target_sensitive=target_sensitive,
        target_extra=target_extra,
        tx=tx,
        ty=ty,
        tz=tz,
        rx=rx,
        ry=ry,
        rz=rz,
        shift_mag_mm=shift_mag_mm,
        rotation_mag_rad=rotation_mag_rad,
    )
    mutableStore = acts.examples.alignment.MutableGeoIdAlignmentStore(geoIdMap)

    # Configure AlignmentDecorator (IOV can be set to global interval)
    cfg = AlignmentDecorator.Config()
    cfg.iovStores = [((0, 1_000_000), mutableStore)]  # Effective for all events
    cfg.garbageCollection = False
    alignDeco = AlignmentDecorator(cfg, acts.logging.INFO)
    s.addContextDecorator(alignDeco)

    logger.info("Reading simulation results from CSV files in {}", str(inputDir))
    s.addReader(
        acts.examples.CsvSimHitReader(
            level=acts.logging.INFO,
            inputDir=str(inputDir),
            inputStem="hits",
            outputSimHits="simhits",
        )
    )

    # Read measurements from CSV
    s.addReader(
        acts.examples.CsvMeasurementReader(
            level=acts.logging.INFO,
            inputDir=str(inputDir),
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="measurement_simhits_map",
            outputMeasurementParticlesMap="measurement_particles_map",
            outputParticleMeasurementsMap="particle_measurements_map",
            inputSimHits="simhits",
        )
    )

    # Read simulated particles from CSV
    s.addReader(
        acts.examples.CsvParticleReader(
            level=acts.logging.INFO,
            inputDir=str(inputDir),
            inputStem="particles_simulated",
            outputParticles="particles_simulated",
        )
    )

    s.addWhiteboardAlias("particles", "particles_simulated")
    s.addWhiteboardAlias("particles_simulated_selected", "particles_simulated")
    s.addWhiteboardAlias("particles_selected", "particles_simulated")

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(2 * u.GeV, None),
            measurements=(5, None),
            removeNeutral=True,
            removeSecondaries=True,
        ),
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_selected",
        initialVarInflation=initialvarinflation,
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        particleHypothesis=acts.ParticleHypothesis.muon,
    )

    # Alignment DOF masks (should match misalignment DOF; local x/y/z)
    Center0 = 1 << 0  # local translation along module local x - matches tx
    Center1 = 1 << 1  # local translation along module local y - matches ty
    Center2 = 1 << 2  # local translation along module local z - matches tz
    Rotation0 = 1 << 3  # rx (rotation around X) - matches rx
    Rotation1 = 1 << 4  # ry (rotation around Y) - matches ry
    Rotation2 = 1 << 5  # rz (rotation around Z) - matches rz

    # Enable alignment DOF based on misalignment DOF
    alignment_dof = 0
    if tx:
        alignment_dof |= Center0
    if ty:
        alignment_dof |= Center1
    if tz:
        alignment_dof |= Center2
    if rx:
        alignment_dof |= Rotation0
    if ry:
        alignment_dof |= Rotation1
    if rz:
        alignment_dof |= Rotation2

    det_elements = alignment_placements
    if len(det_elements) != len(geoIdMap):
        print(
            f"WARNING: misaligned ({len(geoIdMap)}) vs alignment ({len(det_elements)}) count mismatch"
        )

    # Configure AlignmentAlgorithm
    aal_cfg = acts.examples.alignment.AlignmentAlgorithmConfig()
    aal_cfg.inputMeasurements = "measurements"
    aal_cfg.inputProtoTracks = "truth_particle_tracks"
    aal_cfg.inputInitialTrackParameters = "estimatedparameters"
    aal_cfg.outputAlignmentParameters = "alignmentParameters"

    # Set trackingGeometry (for surfaceAccessor)
    aal_cfg.trackingGeometry = trackingGeometry
    aal_cfg.chi2ONdfCutOff = 0.1
    aal_cfg.deltaChi2ONdfCutOff = (5, 0.0001)

    # Create AlignedTransformUpdater, it will update mutableStore after each iteration
    aal_cfg.alignedTransformUpdater = (
        acts.examples.alignment.makeAlignedTransformUpdater(mutableStore)
    )

    aal_cfg.alignedDetElements = det_elements
    aal_cfg.align = acts.examples.alignment.makeAlignmentFunction(
        trackingGeometry, field, acts.logging.INFO
    )
    aal_cfg.maxNumIterations = 100  # Reasonable number of iterations
    aal_cfg.maxNumTracks = 100000
    aal_cfg.iterationState = {
        i: int(alignment_dof) for i in range(aal_cfg.maxNumIterations)
    }

    alignment_algo = acts.examples.alignment.AlignmentAlgorithm(
        aal_cfg, acts.logging.INFO
    )
    s.addAlgorithm(alignment_algo)

    # Run the sequencer
    s.run()

    transform_map = mutableStore.getTransformMap()
    aligned_record = [
        transform_to_alignment_record(gid, trf) for gid, trf in transform_map.items()
    ]

    aligned_file = outputDir / "aligned_transforms.json"
    with open(aligned_file, "w") as f:
        json.dump(aligned_record, f, indent=2)
    print(f"Alignment done: {len(aligned_record)} elements -> {aligned_file.name}")

    return s


# ============================================================================
# Reconstruction functions
# ============================================================================


def runReconstruction(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    inputDir: Path,
    outputDir: Path,
    detector=None,
    numEvents=10,
    reverseFilteringMomThreshold=0 * u.GeV,
    reverseFilteringCovarianceScaling=1,
    s: acts.examples.Sequencer = None,
    mode: str = "nominal",  # "nominal", "misaligned", or "aligned"
    alignmentDir: Path = None,  # Directory containing alignment output files
):
    """
    Run reconstruction with different geometry modes.

    Parameters
    ----------
    mode : str
        - "nominal": Use nominal geometry (no alignment decorator)
        - "misaligned": Use misaligned geometry (loads from misalignment_applied.json)
        - "aligned": Use aligned geometry (loads from aligned_transforms.json)
    alignmentDir : Path
        Directory containing alignment output files (default: alignment_output)
    """
    from acts.examples.simulation import (
        addFatras,
        addDigitization,
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        addKalmanTracks,
    )

    s = s or acts.examples.Sequencer(
        events=numEvents, numThreads=1, logLevel=acts.logging.INFO
    )

    rnd = acts.examples.RandomNumbers(seed=42)
    inputDir = Path(inputDir)
    outputDir = Path(outputDir)

    logger = acts.getDefaultLogger("Reconstruction", acts.logging.INFO)

    # Determine alignment directory
    if alignmentDir is None:
        alignmentDir = Path(__file__).resolve().parent / "alignment_output"

    # Load transforms based on mode
    if mode == "nominal":
        use_alignment_decorator = False
        json_path = None
    elif mode == "misaligned":
        json_path = alignmentDir / "misalignment_applied.json"
        use_alignment_decorator = True
    elif mode == "aligned":
        json_path = alignmentDir / "aligned_transforms.json"
        use_alignment_decorator = True
    else:
        raise ValueError(
            f"Unknown mode: {mode}. Must be 'nominal', 'misaligned', or 'aligned'"
        )

    if use_alignment_decorator:
        if not json_path.exists():
            print(f"Error: Transform file {json_path} does not exist!")
            print(f"Please run alignment mode first to generate the required files.")
            sys.exit(1)

        with open(json_path, "r") as f:
            transform_records = json.load(f)

        geoIdMap = geo_id_map_from_transform_records(transform_records)
        print(
            f"Reconstruction ({mode}): {len(geoIdMap)} aligned elements from {json_path.name}"
        )
        mutableStore = acts.examples.alignment.MutableGeoIdAlignmentStore(geoIdMap)

        cfg = AlignmentDecorator.Config()
        cfg.iovStores = [((0, 1_000_000), mutableStore)]  # Effective for all events
        cfg.garbageCollection = False
        alignDeco = AlignmentDecorator(cfg, acts.logging.INFO)

        s.addContextDecorator(alignDeco)

    logger.info("Reading particles from {}", str(inputDir / "particles.root"))

    s.addReader(
        acts.examples.CsvSimHitReader(
            level=acts.logging.INFO,
            inputDir=str(inputDir),
            inputStem="hits",
            outputSimHits="simhits",
        )
    )

    s.addReader(
        acts.examples.CsvMeasurementReader(
            level=acts.logging.INFO,
            inputDir=str(inputDir),
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="measurement_simhits_map",
            outputMeasurementParticlesMap="measurement_particles_map",
            outputParticleMeasurementsMap="particle_measurements_map",
            inputSimHits="simhits",
        )
    )

    s.addReader(
        acts.examples.CsvParticleReader(
            level=acts.logging.INFO,
            inputDir=str(inputDir),
            inputStem="particles_simulated",
            outputParticles="particles_simulated",
        )
    )

    s.addWhiteboardAlias("particles", "particles_simulated")
    s.addWhiteboardAlias("particles_simulated_selected", "particles_simulated")
    s.addWhiteboardAlias("particles_selected", "particles_simulated")

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(2 * u.GeV, None),
            measurements=(5, None),
            removeNeutral=True,
            removeSecondaries=True,
        ),
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_selected",
        initialVarInflation=[1] * 6,
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        particleHypothesis=acts.ParticleHypothesis.muon,
    )

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
        reverseFilteringMomThreshold,
        reverseFilteringCovarianceScaling,
    )

    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected-tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=5,
            ),
        )
    )
    s.addWhiteboardAlias("tracks", "selected-tracks")

    s.addWriter(
        acts.examples.root.RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_kf.root"),
        )
    )

    s.addWriter(
        acts.examples.root.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf.root"),
        )
    )

    s.addWriter(
        acts.examples.root.RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_kf.root"),
        )
    )

    return s


# ============================================================================
# Main entry point
# ============================================================================

if "__main__" == __name__:
    parser = argparse.ArgumentParser(
        description="Complete simulation, alignment, and reconstruction workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
            # Run complete workflow (simulation -> alignment -> reconstruction with all geometries)
            python alignment_odd.py full

            # Run individual steps
            python alignment_odd.py simulation
            python alignment_odd.py alignment
            python alignment_odd.py reconstruction --mode nominal
            python alignment_odd.py reconstruction --mode misaligned
            python alignment_odd.py reconstruction --mode aligned
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to execute")

    # Full workflow command
    full_parser = subparsers.add_parser(
        "full",
        help="Run complete workflow: simulation -> alignment -> reconstruction (all geometries)",
    )
    full_parser.add_argument(
        "--num-events", type=int, default=100000, help="Number of events"
    )
    full_parser.add_argument(
        "--simulation-dir",
        type=str,
        default=None,
        help="Simulation output directory (default: simulation_output)",
    )
    full_parser.add_argument(
        "--alignment-dir",
        type=str,
        default=None,
        help="Alignment output directory (default: alignment_output)",
    )
    full_parser.add_argument(
        "--reconstruction-dir",
        type=str,
        default=None,
        help="Reconstruction output directory (default: reconstruction_output)",
    )

    # Simulation command
    sim_parser = subparsers.add_parser("simulation", help="Run simulation")
    sim_parser.add_argument(
        "--num-events", type=int, default=100000, help="Number of events"
    )
    sim_parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: simulation_output)",
    )

    # Alignment command
    align_parser = subparsers.add_parser("alignment", help="Run alignment algorithm")
    align_parser.add_argument(
        "--num-events", type=int, default=100000, help="Number of events"
    )
    align_parser.add_argument(
        "--input-dir",
        type=str,
        default=None,
        help="Input directory (default: simulation_output)",
    )
    align_parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: alignment_output)",
    )

    # Reconstruction command
    recon_parser = subparsers.add_parser("reconstruction", help="Run reconstruction")
    recon_parser.add_argument(
        "--mode",
        type=str,
        choices=["nominal", "misaligned", "aligned"],
        default="nominal",
        help="Reconstruction mode",
    )
    recon_parser.add_argument(
        "--num-events", type=int, default=100000, help="Number of events"
    )
    recon_parser.add_argument(
        "--input-dir",
        type=str,
        default=None,
        help="Input directory (default: simulation_output)",
    )
    recon_parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: reconstruction_output)",
    )

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # OpenDataDetector
    from acts.examples.odd import getOpenDataDetector

    def build_odd(misaligned: bool = False):
        detector = getOpenDataDetector(misaligned=misaligned)
        return detector, detector.trackingGeometry()

    # Setup common paths and detector
    digiConfigFile = srcdir / "Examples/Configs/odd-digi-smearing-config.json"
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    # Execute commands
    if args.command == "full":
        # Full workflow: simulation -> alignment -> reconstruction (all three geometries)
        print("\n" + "=" * 80)
        print("=" * 80)
        print("RUNNING COMPLETE WORKFLOW")
        print("=" * 80)
        print("=" * 80)

        simulationDir = (
            Path(args.simulation_dir)
            if args.simulation_dir
            else Path.cwd() / "simulation_output"
        )
        alignmentDir = (
            Path(args.alignment_dir)
            if args.alignment_dir
            else Path.cwd() / "alignment_output"
        )
        reconstructionDir = (
            Path(args.reconstruction_dir)
            if args.reconstruction_dir
            else Path.cwd() / "reconstruction_output"
        )

        # Step 1: Simulation (using nominal geometry)
        print("\n" + "=" * 80)
        print("STEP 1: SIMULATION (using NOMINAL geometry)")
        print("=" * 80)
        detector_nom, trackingGeometry_nom = build_odd(misaligned=False)

        runSimulation(
            trackingGeometry=trackingGeometry_nom,
            field=field,
            digiConfigFile=digiConfigFile,
            outputDir=simulationDir,
            detector=detector_nom,
            numEvents=args.num_events,
        ).run()

        # Step 2: Alignment (using misaligned geometry)
        print("\n" + "=" * 80)
        print("STEP 2: ALIGNMENT (using MISALIGNED geometry)")
        print("=" * 80)
        alignmentDir.mkdir(exist_ok=True)

        detector_mis, trackingGeometry_mis = build_odd(misaligned=True)

        runAlignment(
            trackingGeometry=trackingGeometry_mis,
            field=field,
            digiConfigFile=digiConfigFile,
            inputDir=simulationDir,
            outputDir=alignmentDir,
            detector=detector_mis,
            numEvents=args.num_events,
        )

        # Step 3: Reconstruction with nominal geometry
        print("\n" + "=" * 80)
        print("STEP 3: RECONSTRUCTION WITH NOMINAL GEOMETRY")
        print("=" * 80)
        # detector_nom and trackingGeometry_nom already created in Step 1
        reconDir_nom = reconstructionDir / "nominal"
        reconDir_nom.mkdir(parents=True, exist_ok=True)

        runReconstruction(
            trackingGeometry=trackingGeometry_nom,
            field=field,
            digiConfigFile=digiConfigFile,
            inputDir=simulationDir,
            outputDir=reconDir_nom,
            detector=detector_nom,
            numEvents=args.num_events,
            mode="nominal",
            alignmentDir=alignmentDir,
        ).run()

        # Step 4: Reconstruction with misaligned geometry
        print("\n" + "=" * 80)
        print("STEP 4: RECONSTRUCTION WITH MISALIGNED GEOMETRY")
        print("=" * 80)
        reconDir_mis = reconstructionDir / "misaligned"
        reconDir_mis.mkdir(parents=True, exist_ok=True)

        runReconstruction(
            trackingGeometry=trackingGeometry_mis,
            field=field,
            digiConfigFile=digiConfigFile,
            inputDir=simulationDir,
            outputDir=reconDir_mis,
            detector=detector_mis,
            numEvents=args.num_events,
            mode="misaligned",
            alignmentDir=alignmentDir,
        ).run()

        # Step 5: Reconstruction with aligned geometry
        print("\n" + "=" * 80)
        print("STEP 5: RECONSTRUCTION WITH ALIGNED GEOMETRY")
        print("=" * 80)
        reconDir_ali = reconstructionDir / "aligned"
        reconDir_ali.mkdir(parents=True, exist_ok=True)

        runReconstruction(
            trackingGeometry=trackingGeometry_mis,
            field=field,
            digiConfigFile=digiConfigFile,
            inputDir=simulationDir,
            outputDir=reconDir_ali,
            detector=detector_mis,
            numEvents=args.num_events,
            mode="aligned",
            alignmentDir=alignmentDir,
        ).run()

        print("\n" + "=" * 80)
        print("=" * 80)
        print("COMPLETE WORKFLOW FINISHED SUCCESSFULLY!")
        print("=" * 80)
        print("=" * 80)
        print(f"\nOutput directories:")
        print(f"  - Simulation: {simulationDir}")
        print(f"  - Alignment: {alignmentDir}")
        print(f"  - Reconstruction (nominal): {reconDir_nom}")
        print(f"  - Reconstruction (misaligned): {reconDir_mis}")
        print(f"  - Reconstruction (aligned): {reconDir_ali}")

    elif args.command == "simulation":
        outputDir = (
            Path(args.output_dir)
            if args.output_dir
            else Path.cwd() / "simulation_output"
        )
        detector, trackingGeometry = build_odd(misaligned=True)

        runSimulation(
            trackingGeometry=trackingGeometry,
            field=field,
            digiConfigFile=digiConfigFile,
            outputDir=outputDir,
            detector=detector,
            numEvents=args.num_events,
        ).run()

    elif args.command == "alignment":
        inputDir = (
            Path(args.input_dir) if args.input_dir else Path.cwd() / "simulation_output"
        )
        outputDir = (
            Path(args.output_dir)
            if args.output_dir
            else Path.cwd() / "alignment_output"
        )
        outputDir.mkdir(exist_ok=True)

        if not inputDir.exists():
            print(f"Error: Input directory {inputDir} does not exist!")
            print("Please run simulation first to generate simulation data.")
            sys.exit(1)

        detector, trackingGeometry = build_odd(misaligned=True)

        runAlignment(
            trackingGeometry=trackingGeometry,
            field=field,
            digiConfigFile=digiConfigFile,
            inputDir=inputDir,
            outputDir=outputDir,
            detector=detector,
            numEvents=args.num_events,
        )

    elif args.command == "reconstruction":
        inputDir = (
            Path(args.input_dir) if args.input_dir else Path.cwd() / "simulation_output"
        )
        outputDir = (
            Path(args.output_dir)
            if args.output_dir
            else Path.cwd() / "reconstruction_output"
        )
        outputDir.mkdir(exist_ok=True)

        if not inputDir.exists():
            print(f"Error: Input directory {inputDir} does not exist!")
            print("Please run simulation first to generate simulation data.")
            sys.exit(1)

        # Choose detector based on mode
        misaligned = args.mode != "nominal"
        detector, trackingGeometry = build_odd(misaligned=misaligned)

        runReconstruction(
            trackingGeometry=trackingGeometry,
            field=field,
            digiConfigFile=digiConfigFile,
            inputDir=inputDir,
            outputDir=outputDir,
            detector=detector,
            numEvents=args.num_events,
            mode=args.mode,
        ).run()

        print("\n" + "=" * 60)
        print("Reconstruction completed successfully!")
        print("=" * 60)
        print(f"\nOutput files saved to: {outputDir}")
        print("\nGenerated files:")
        print("  - trackstates_kf.root: Track states")
        print("  - tracksummary_kf.root: Track summary")
        print("  - performance_kf.root: Reconstruction performance")
