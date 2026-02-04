from pathlib import Path
from typing import Optional
import json
import random
import numpy as np
import ctypes
import os
import argparse
import sys

import acts
import acts.examples
import acts.examples.alignment
from acts.examples.alignment import AlignmentDecorator

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

    logger = acts.logging.getLogger("Simulation")

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
    logger.info("Simulation output will be saved to: %s", outputDir)
    logger.info(
        "particle_measurements_map will be generated after simulation completes"
    )

    return s


# ============================================================================
# Alignment functions
# ============================================================================


def setupMisalignment(
    detector_elements: list,
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
    Setup misalignment for detector elements

    Parameters
    ----------
    detector_elements : list
        List of all detector elements from odd_transforms.json
    outputDir : Path
        Directory to save misalignment record
    target_volume : int or list of ints
        Filter for volume IDs (-1 means no restriction, or list of specific IDs)
    target_layer : int or list of ints
        Filter for layer IDs (-1 means no restriction, or list of specific IDs)
    target_sensitive : int or list of ints
        Filter for sensitive IDs (-1 means no restriction, or list of specific IDs)
    target_extra : int
        Filter for extra ID (-1 means no restriction)
    tx, ty, tz : int
        Translation DOF flags (0 or 1)
    rx, ry, rz : int
        Rotation DOF flags (0 or 1)
    shift_mag_mm : float
        Random shift magnitude in mm
    rotation_mag_rad : float
        Random rotation magnitude in radians

    Returns
    -------
    tuple
        (geoIdMap, selected_elements, misalignment_record)
    """
    print("\n" + "=" * 60)
    print("Setting up misalignment in LOCAL coordinate system...")
    print("=" * 60)

    print(f"Misalignment DOF enabled (in local coordinates):")
    print(f"  Translation: tx={tx} (local X), ty={ty} (local Y), tz={tz} (local Z)")
    print(f"  Rotation:    rx={rx} (local X), ry={ry} (local Y), rz={rz} (local Z)")
    print(
        f"  Magnitudes:  ±{shift_mag_mm} mm (translation), ±{rotation_mag_rad} rad (rotation)"
    )
    print(f"\nNote: Using AlignmentGenerator LocalShift and LocalRotation:")
    print(f"  - LocalShift: Applies translation along specified local axis")
    print(f"  - LocalRotation: Applies rotation around specified local axis")
    print(f"  - Transformations are automatically handled in local coordinate system")

    # Filter elements to misalign
    selected_elements = []
    for element in detector_elements:
        geoid_parts = parse_geoid(element["ID"])
        if match_filter(
            geoid_parts, target_volume, target_layer, target_sensitive, target_extra
        ):
            selected_elements.append(element)

    print(f"\nFilter criteria:")

    # Format volume output (handle hierarchical dict, list, or single value)
    if isinstance(target_volume, dict):
        volume_str = "Hierarchical config:\n"
        for vol_id, layer_config in sorted(target_volume.items()):
            if layer_config == -1:
                volume_str += f"    Volume {vol_id}: ALL layers, ALL sensors\n"
            elif isinstance(layer_config, dict):
                volume_str += f"    Volume {vol_id}:\n"
                for lay_id, sensor_config in sorted(layer_config.items()):
                    if sensor_config == -1:
                        volume_str += f"      Layer {lay_id}: ALL sensors\n"
                    elif isinstance(sensor_config, list):
                        sensor_list = sensor_config[:10]
                        sensor_str = (
                            f"{sensor_list}{'...' if len(sensor_config) > 10 else ''}"
                        )
                        volume_str += f"      Layer {lay_id}: sensors {sensor_str}\n"
                    else:
                        volume_str += f"      Layer {lay_id}: sensor {sensor_config}\n"
            elif isinstance(layer_config, list):
                volume_str += (
                    f"    Volume {vol_id}: layers {layer_config}, ALL sensors\n"
                )
            else:
                volume_str += (
                    f"    Volume {vol_id}: layer {layer_config}, ALL sensors\n"
                )
        volume_str = volume_str.rstrip()
    elif target_volume == -1:
        volume_str = "ALL"
    elif isinstance(target_volume, list):
        volume_str = (
            f"{target_volume}"
            if len(target_volume) <= 5
            else f"{len(target_volume)} IDs: {target_volume[:5]}..."
        )
    else:
        volume_str = str(target_volume)
    print(f"  Volume:    {volume_str}")

    # Format layer output (only used for backward compatibility)
    if isinstance(target_volume, dict):
        layer_str = "N/A (using hierarchical config)"
    elif target_layer == -1:
        layer_str = "ALL"
    elif isinstance(target_layer, list):
        layer_str = (
            f"{target_layer}"
            if len(target_layer) <= 5
            else f"{len(target_layer)} IDs: {target_layer[:5]}..."
        )
    else:
        layer_str = str(target_layer)
    print(f"  Layer:     {layer_str}")

    # Format sensitive output (only used for backward compatibility)
    if isinstance(target_volume, dict):
        sensitive_str = "N/A (using hierarchical config)"
    elif target_sensitive == -1:
        sensitive_str = "ALL"
    elif isinstance(target_sensitive, list):
        sensitive_str = f'{len(target_sensitive)} IDs: {target_sensitive[:10]}{"..." if len(target_sensitive) > 10 else ""}'
    else:
        sensitive_str = str(target_sensitive)
    print(f"  Sensitive: {sensitive_str}")

    print(f"  Extra:     {target_extra if target_extra != -1 else 'ALL'}")
    print(f"Selected elements: {len(selected_elements)}")
    print(f"Shift magnitude: ±{shift_mag_mm} mm")
    print(f"Rotation magnitude: ±{rotation_mag_rad} rad")

    # Construct geoIdMap, apply random misalignment in local coordinate system
    # Also build misalignment record simultaneously
    # Using AlignmentGenerator LocalShift and LocalRotation for proper local transformations
    geoIdMap = {}
    misalignment_record = []

    # Random shift with asymmetric distribution
    # Each axis randomly chooses positive or negative side
    def asymmetric_random(mag, enabled):
        if not enabled:
            return 0.0
        if random.random() < 0.5:
            return random.uniform(-mag * 1.5, -mag * 0.5)
        else:
            return random.uniform(mag * 0.5, mag * 1.5)

    for element in selected_elements:
        geoid_parts = parse_geoid(element["ID"])
        nom_trans = element["Translation"]
        nom_rot_matrix = element["RotationMatrix"]

        # Create GeometryIdentifier
        gid_kwargs = {
            "volume": geoid_parts.get("vol", 0),
            "layer": geoid_parts.get("lay", 0),
            "sensitive": geoid_parts.get("sen", 0),
        }
        if "ext" in geoid_parts:
            gid_kwargs["extra"] = geoid_parts["ext"]
        gid = acts.GeometryIdentifier(**gid_kwargs)

        # Read rotation matrix directly from JSON
        R_nominal = np.array(nom_rot_matrix)

        # Create nominal Transform3
        nom_rotation = acts.RotationMatrix3(
            acts.Vector3(R_nominal[0, 0], R_nominal[1, 0], R_nominal[2, 0]),
            acts.Vector3(R_nominal[0, 1], R_nominal[1, 1], R_nominal[2, 1]),
            acts.Vector3(R_nominal[0, 2], R_nominal[1, 2], R_nominal[2, 2]),
        )
        trf = acts.Transform3(
            acts.Vector3(nom_trans[0], nom_trans[1], nom_trans[2]), nom_rotation
        )

        # Store local shift and rotation values for recording
        shift_local_x = asymmetric_random(shift_mag_mm, tx != 0)
        shift_local_y = asymmetric_random(shift_mag_mm, ty != 0)
        shift_local_z = 0
        rot_local_z = asymmetric_random(rotation_mag_rad, rz != 0)
        rot_local_y = 0
        rot_local_x = 0

        # ===== Apply LocalShift using AlignmentGenerator =====
        # LocalShift applies translation in local coordinate system
        if tx != 0:
            lShiftX = acts.examples.alignment.AlignmentGeneratorLocalShift()
            lShiftX.axisDirection = acts.AxisDirection.AxisX
            lShiftX.shift = shift_local_x
            lShiftX(trf)  # Apply local shift along X axis

        if ty != 0:
            lShiftY = acts.examples.alignment.AlignmentGeneratorLocalShift()
            lShiftY.axisDirection = acts.AxisDirection.AxisY
            lShiftY.shift = shift_local_y
            lShiftY(trf)  # Apply local shift along Y axis

        # ===== Apply LocalRotation using AlignmentGenerator =====
        # LocalRotation applies rotation around local axes
        # IMPORTANT: The axis parameter must be the local axis direction in GLOBAL coordinates
        # We need to extract the local axes from the current transform's rotation matrix
        current_rotation = trf.rotation

        if rz != 0:
            lRotZ = acts.examples.alignment.AlignmentGeneratorLocalRotation()
            # Local Z axis in global coordinates = third column of rotation matrix
            local_z_axis = acts.Vector3(
                current_rotation[0, 2], current_rotation[1, 2], current_rotation[2, 2]
            )
            lRotZ.axis = local_z_axis
            lRotZ.angle = rot_local_z
            lRotZ(trf)  # Apply rotation around local Z axis
            # Update rotation matrix after applying rotation
            current_rotation = trf.rotation

        # Store the misaligned transform
        geoIdMap[gid] = trf

        # Extract misaligned transform components
        misaligned_trans = trf.translation
        misaligned_rot = trf.rotation

        # Extract rotation matrix as nested list
        rotation_matrix = [
            [
                float(misaligned_rot[0, 0]),
                float(misaligned_rot[0, 1]),
                float(misaligned_rot[0, 2]),
            ],
            [
                float(misaligned_rot[1, 0]),
                float(misaligned_rot[1, 1]),
                float(misaligned_rot[1, 2]),
            ],
            [
                float(misaligned_rot[2, 0]),
                float(misaligned_rot[2, 1]),
                float(misaligned_rot[2, 2]),
            ],
        ]

        # Simplified record: only ID, Translation, and RotationMatrix
        record = {
            "ID": element["ID"],
            "Translation": [
                misaligned_trans[0],
                misaligned_trans[1],
                misaligned_trans[2],
            ],
            "RotationMatrix": rotation_matrix,
        }
        misalignment_record.append(record)

    print(f"\nTotal misaligned detector elements: {len(geoIdMap)}")

    # Save misalignment record
    misalignment_file = outputDir / "misalignment_applied.json"
    with open(misalignment_file, "w") as f:
        json.dump(misalignment_record, f, indent=2)
    print(f"Misalignment record saved to: {misalignment_file}")

    return geoIdMap, selected_elements, misalignment_record


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

    logger = acts.logging.getLogger("Alignment")

    # Load odd_transforms.json
    json_path = Path(__file__).resolve().parent / "Configs/odd_transforms.json"
    print(f"Loading detector transforms from: {json_path}")

    with open(json_path, "r") as f:
        detector_elements = json.load(f)
    print(f"Total detector elements in JSON: {len(detector_elements)}")

    # Setup misalignment parameters
    # Format 1: Hierarchical structure - specify volume -> layer -> sensors
    # Example: {16: {4: [1, 2, 3], 5: [1, 2]}, 17: {2: [1]}}
    #   - Volume 16, Layer 4: sensors [1, 2, 3]
    #   - Volume 16, Layer 5: sensors [1, 2]
    #   - Volume 17, Layer 2: sensor [1]
    # Use -1 for layer value to match all sensors in that layer
    # Use -1 for sensor list to match all sensors: {16: {4: -1}} means all sensors in vol16/lay4
    # Format 1: Use hierarchical structure (dict) - specify volume -> layer -> sensors
    target_volume = {24: {4: -1}}
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

    geoIdMap, selected_elements, misalignment_record = setupMisalignment(
        detector_elements=detector_elements,
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

    logger.info("Reading simulation results from CSV files in %s", inputDir)
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

    print("\n" + "=" * 60)
    print("Setting up alignment algorithm...")
    print("=" * 60)

    # Alignment DOF masks (should match misalignment DOF)
    Center0 = 1 << 0  # dx (translation X) - matches tx
    Center1 = 1 << 1  # dy (translation Y) - matches ty
    Center2 = 1 << 2  # dz (translation Z) - matches tz
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

    print(f"Alignment DOF enabled (matches misalignment):")
    print(
        f"  dx={bool(alignment_dof & Center0)}, dy={bool(alignment_dof & Center1)}, dz={bool(alignment_dof & Center2)}"
    )
    print(
        f"  rx={bool(alignment_dof & Rotation0)}, ry={bool(alignment_dof & Rotation1)}, rz={bool(alignment_dof & Rotation2)}"
    )
    print(f"  Alignment mask value: {alignment_dof} (binary: {bin(alignment_dof)})")

    # Collect detector elements to align (same as misaligned elements)
    # Use the same filtering logic as misalignment section
    print("\nCollecting detector elements for alignment...")
    print(f"Using same filter as misalignment:")
    if isinstance(target_volume, dict):
        print(f"  Volume:    Hierarchical config (see misalignment section above)")
    else:
        print(f"  Volume:    {target_volume if target_volume != -1 else 'ALL'}")
    if not isinstance(target_volume, dict):
        print(f"  Layer:     {target_layer if target_layer != -1 else 'ALL'}")
        print(f"  Sensitive: {target_sensitive if target_sensitive != -1 else 'ALL'}")
    print(f"  Extra:     {target_extra if target_extra != -1 else 'ALL'}")

    det_elements = []

    def visit(surface: acts.Surface) -> bool:
        # Use surfacePlacement instead of associatedDetectorElement
        # because DetectorElementBase is not registered in Python bindings
        placement = acts.examples.alignment.surfacePlacement(surface)
        if placement is not None:
            gid = surface.geometryId
            # Use same filter as misalignment
            if match_surface_to_filter(
                gid, target_volume, target_layer, target_sensitive, target_extra
            ):
                det_elements.append(placement)
        return True

    trackingGeometry.visitSurfaces(visit)

    print(f"Total detector elements selected for alignment: {len(det_elements)}")

    # Verify the count matches misaligned elements
    if len(det_elements) != len(geoIdMap):
        print(
            f"WARNING: Mismatch between misaligned ({len(geoIdMap)}) and alignment ({len(det_elements)}) element counts!"
        )
    else:
        print(f"Alignment element count matches misaligned elements")

    # Configure AlignmentAlgorithm
    aal_cfg = acts.examples.alignment.AlignmentAlgorithmConfig()
    aal_cfg.inputMeasurements = "measurements"
    aal_cfg.inputProtoTracks = "truth_particle_tracks"
    aal_cfg.inputInitialTrackParameters = "estimatedparameters"
    aal_cfg.outputAlignmentParameters = "alignmentParameters"

    # Set trackingGeometry (for surfaceAccessor)
    aal_cfg.trackingGeometry = trackingGeometry
    aal_cfg.chi2ONdfCutOff = 0.1
    aal_cfg.deltaChi2ONdfCutOff = (5, 0.00001)

    # Create AlignedTransformUpdater, it will update mutableStore after each iteration
    aal_cfg.alignedTransformUpdater = acts.examples.alignment.makeAlignedTransformUpdater(
        mutableStore
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

    alignment_algo = acts.examples.alignment.AlignmentAlgorithm(aal_cfg, acts.logging.INFO)
    s.addAlgorithm(alignment_algo)

    # Run the sequencer
    s.run()

    print("\n" + "=" * 60)
    print("Extracting aligned transforms...")
    print("=" * 60)

    # Extract aligned transforms from mutableStore after alignment
    aligned_record = []
    transform_map = mutableStore.getTransformMap()
    gctx = acts.GeometryContext()

    # Collect surfaces that were aligned (same filter as misalignment)
    aligned_surfaces = []

    def visit(surface: acts.Surface) -> bool:
        # Use surfacePlacement to check if surface has a placement
        placement = acts.examples.alignment.surfacePlacement(surface)
        if placement is not None:
            gid = surface.geometryId
            # Use same filter as misalignment
            if match_surface_to_filter(
                gid, target_volume, target_layer, target_sensitive, target_extra
            ):
                aligned_surfaces.append(surface)
        return True

    trackingGeometry.visitSurfaces(visit)

    for surface in aligned_surfaces:
        gid = surface.geometryId
        element_id = geoid_to_id_string(gid)

        # Find the matching GeometryIdentifier in transform_map
        aligned_trf = None
        for map_gid, trf in transform_map.items():
            if str(map_gid) == str(gid):
                aligned_trf = trf
                break

        if aligned_trf is not None:
            aligned_trans = aligned_trf.translation
            aligned_rot = aligned_trf.rotation

            # Extract aligned rotation matrix as nested list
            aligned_rotation_matrix = [
                [
                    float(aligned_rot[0, 0]),
                    float(aligned_rot[0, 1]),
                    float(aligned_rot[0, 2]),
                ],
                [
                    float(aligned_rot[1, 0]),
                    float(aligned_rot[1, 1]),
                    float(aligned_rot[1, 2]),
                ],
                [
                    float(aligned_rot[2, 0]),
                    float(aligned_rot[2, 1]),
                    float(aligned_rot[2, 2]),
                ],
            ]

            # Simplified record: only ID, Translation, and RotationMatrix
            record = {
                "ID": element_id,
                "Translation": [aligned_trans[0], aligned_trans[1], aligned_trans[2]],
                "RotationMatrix": aligned_rotation_matrix,
            }
            aligned_record.append(record)
        else:
            print(f"Warning: Could not find aligned transform for {element_id}")

    # Save aligned transforms record
    aligned_file = outputDir / "aligned_transforms.json"
    with open(aligned_file, "w") as f:
        json.dump(aligned_record, f, indent=2)
    print(f"\nAligned transforms saved to: {aligned_file}")
    print(f"Total aligned elements: {len(aligned_record)}")

    # Summary information
    print(f"\nTotal aligned elements: {len(aligned_record)}")

    print("\n" + "=" * 60)
    print("Alignment completed successfully!")
    print("=" * 60)
    print(f"\nOutput files saved to: {outputDir}")
    print("\nGenerated files:")
    print("  - misalignment_applied.json: Applied misalignment transforms")
    print("    * Includes Translation and RotationMatrix")
    print("  - aligned_transforms.json: Final aligned transforms")
    print("    * Includes Translation and RotationMatrix")

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

    logger = acts.logging.getLogger("Reconstruction")

    # Determine alignment directory
    if alignmentDir is None:
        alignmentDir = Path(__file__).resolve().parent / "alignment_output"

    # Load transforms based on mode
    if mode == "nominal":
        print("=" * 60)
        print("Running reconstruction with NOMINAL geometry...")
        print("=" * 60)
        # No alignment decorator needed for nominal geometry
        use_alignment_decorator = False
    elif mode == "misaligned":
        print("=" * 60)
        print("Running reconstruction with MISALIGNED geometry...")
        print("=" * 60)
        json_path = alignmentDir / "misalignment_applied.json"
        print(f"Loading misaligned transforms from: {json_path}")
        use_alignment_decorator = True
    elif mode == "aligned":
        print("=" * 60)
        print("Running reconstruction with ALIGNED geometry...")
        print("=" * 60)
        json_path = alignmentDir / "aligned_transforms.json"
        print(f"Loading aligned transforms from: {json_path}")
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
            detector_elements = json.load(f)

        print(f"Total detector elements in JSON: {len(detector_elements)}")
        print(f"Using all elements: {len(detector_elements)}")
        print(f"Loading transforms from JSON file...")
        print()

        geoIdMap = {}
        for idx, element in enumerate(detector_elements):
            geoid_parts = parse_geoid(element["ID"])

            gid_kwargs = {
                "volume": geoid_parts.get("vol", 0),
                "layer": geoid_parts.get("lay", 0),
                "sensitive": geoid_parts.get("sen", 0),
            }
            if "ext" in geoid_parts:
                gid_kwargs["extra"] = geoid_parts["ext"]
            gid = acts.GeometryIdentifier(**gid_kwargs)

            # Use transform data directly from JSON (simplified structure)
            trans_data = element["Translation"]
            rot_matrix = element["RotationMatrix"]

            trans = acts.Vector3(trans_data[0], trans_data[1], trans_data[2])

            R = np.array(rot_matrix)
            rotation = acts.RotationMatrix3(
                acts.Vector3(R[0, 0], R[1, 0], R[2, 0]),
                acts.Vector3(R[0, 1], R[1, 1], R[2, 1]),
                acts.Vector3(R[0, 2], R[1, 2], R[2, 2]),
            )

            trf = acts.Transform3(trans, rotation)
            geoIdMap[gid] = trf

        print(f"\nTotal detector elements loaded: {len(geoIdMap)}")
        mutableStore = acts.examples.alignment.MutableGeoIdAlignmentStore(geoIdMap)

        cfg = AlignmentDecorator.Config()
        cfg.iovStores = [((0, 1_000_000), mutableStore)]  # Effective for all events
        cfg.garbageCollection = False
        alignDeco = AlignmentDecorator(cfg, acts.logging.INFO)

        s.addContextDecorator(alignDeco)
    else:
        print("Using nominal geometry (no alignment decorator)")

    logger.info("Reading particles from %s", inputDir / "particles.root")

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
        acts.examples.RootTrackStatesWriter(
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
        acts.examples.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf.root"),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
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
        detector_nom = getOpenDataDetector(misaligned=False)
        trackingGeometry_nom = detector_nom.trackingGeometry()

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

        detector_mis = getOpenDataDetector(misaligned=True)
        trackingGeometry_mis = detector_mis.trackingGeometry()

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
        detector = getOpenDataDetector(misaligned=True)
        trackingGeometry = detector.trackingGeometry()

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

        detector = getOpenDataDetector(misaligned=True)
        trackingGeometry = detector.trackingGeometry()

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
        if args.mode == "nominal":
            detector = getOpenDataDetector(misaligned=False)
        else:
            detector = getOpenDataDetector(misaligned=True)

        trackingGeometry = detector.trackingGeometry()

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
