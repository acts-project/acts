#!/usr/bin/env python3
"""
ODD full chain with the track reconstruction running on a device through traccc

    simulation -> digitization (ACTS, host)
    -> measurements to traccc -> host to device copy -> sorting
    -> spacepoint formation -> (triplet) seeding -> parameter estimation
    -> CKF
    -> device to host copy -> tracks to ACTS
    -> truth matching + performance writers (ACTS, host)

The ODD (Gen3) tracking geometry is converted to a detray detector in memory,
no detray JSON files are needed. The device algorithms are backend agnostic,
the backend is picked with --backend.
"""

import argparse
import json
from pathlib import Path

import acts
import acts.examples
import acts.covfie
import acts.detray
import acts.json
import acts.vecmem
import acts.examples.traccc as traccc
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory
from acts.examples.simulation import (
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    addParticleGun,
    addPythia8,
    addGenParticleSelection,
    addFatras,
    addSimParticleSelection,
    addDigitization,
    addDigiParticleSelection,
)

u = acts.UnitConstants

actsDir = Path(__file__).parent.parent.parent.parent


def classifyVolume(name):
    """Classify a Gen3 ODD sensitive volume by name into pixel, short strip
    or long strip."""
    lower = name.lower()
    if "pixel" in lower:
        return "pixel"
    if "shortstrip" in lower or "sstrip" in lower:
        return "sstrip"
    if "longstrip" in lower or "lstrip" in lower:
        return "lstrip"
    return None


def makeGen3SmearingConfig(trackingGeometry, outputFile):
    """Write a Gen3 version of the ODD smearing digitization config

    The Gen1 config is keyed by Gen1 volume ids. The Gen3 sensitive volumes
    are classified by name, pixels get the 2D Gen1 pixel smearing, strips
    only the precise (loc0) coordinate. Only pixel measurements are 2D and
    thus turned into spacepoints by traccc.
    """
    gen1 = json.loads(
        (actsDir / "Examples/Configs/odd-digi-smearing-config.json").read_text()
    )
    params = {e["volume"]: e["value"] for e in gen1["entries"]}

    def localOnly(value, indices):
        return {"smearing": [s for s in value["smearing"] if s["index"] in indices]}

    smearing = {
        "pixel": localOnly(params[17], (0, 1)),
        "sstrip": localOnly(params[24], (0,)),
        "lstrip": localOnly(params[29], (0,)),
    }

    volumes = {}

    def visit(surface):
        geoId = surface.geometryId
        if geoId.sensitive == 0:
            return
        volumes.setdefault(geoId.volume, None)

    trackingGeometry.visitSurfaces(visit)

    class VolumeNameVisitor(acts.TrackingGeometryMutableVisitor):
        def __init__(self):
            super().__init__()
            self.names = {}

        def visitVolume(self, volume):
            self.names[volume.geometryId.volume] = volume.volumeName

    visitor = VolumeNameVisitor()
    trackingGeometry.apply(visitor)
    names = visitor.names

    entries = []
    for volumeId in sorted(volumes):
        kind = classifyVolume(names.get(volumeId, ""))
        if kind is None:
            raise RuntimeError(
                f"Cannot classify sensitive volume {volumeId} "
                f"'{names.get(volumeId)}' for the digitization config"
            )
        entries.append({"volume": volumeId, "value": smearing[kind]})

    config = {
        "acts-geometry-hierarchy-map": {
            "format-version": 0,
            "value-identifier": "digitization-configuration",
        },
        "entries": entries,
    }
    outputFile.write_text(json.dumps(config, indent=2))
    return outputFile


def makeDetrayDetector(trackingGeometry, gContext, logLevel):
    """Convert the ACTS tracking geometry into a detray detector in memory,
    together with the ACTS <-> detray surface identifier map"""
    mr = acts.vecmem.HostMemoryResource()

    payloadConfig = acts.detray.DetrayPayloadConverter.Config()
    payloadConfig.beampipeVolume = trackingGeometry.findVolumeByName("BeamPipe")
    payloadConfig.convertMaterial = True
    payloadConfig.convertSurfaceGrids = True
    payloadConverter = acts.detray.DetrayPayloadConverter(payloadConfig, logLevel)

    converterConfig = acts.detray.DetrayGeometryConverter.Config()
    converterConfig.payloadConverter = payloadConverter
    converter = acts.detray.DetrayGeometryConverter(converterConfig, logLevel)

    # The detector keeps the memory resource alive
    return traccc.convertToTracccOddDetector(converter, mr, gContext, trackingGeometry)


def makeBackend(name, logLevel):
    if name == "cuda":
        if not hasattr(traccc, "CudaTracccBackend"):
            raise RuntimeError("ACTS was built without the traccc CUDA backend")
        return traccc.CudaTracccBackend(device=0, level=logLevel)
    raise ValueError(f"Unknown traccc backend '{name}'")


def addTracccDeviceReconstruction(
    s, detector, magneticField, copyTool, logLevel, bFieldInZ
):
    """Add the traccc device chain from ACTS measurements to ACTS tracks"""
    s.addAlgorithm(
        traccc.ActsMeasToTracccAlg(
            level=logLevel,
            detector=detector,
            inputMeasurements="measurements",
            outputMeasurements="traccc_host_measurements",
        )
    )

    s.addAlgorithm(
        traccc.TracccHostToDeviceMeasurementCopyAlg(
            level=logLevel,
            copyTool=copyTool,
            input="traccc_host_measurements",
            output="traccc_device_measurements_unsorted",
        )
    )

    s.addAlgorithm(
        traccc.TracccDeviceMeasurementSortingAlg(
            level=logLevel,
            detector=detector,
            inputMeasurements="traccc_device_measurements_unsorted",
            outputMeasurements="traccc_device_measurements",
        )
    )

    s.addAlgorithm(
        traccc.TracccDeviceSpacepointFormationAlg(
            level=logLevel,
            detector=detector,
            inputMeasurements="traccc_device_measurements",
            outputSpacepoints="traccc_device_spacepoints",
        )
    )

    # example config settings for seed finder
    seedFinder = traccc.SeedFinderConfig()
    seedFinder.collisionRegionMin = -150 * u.mm
    seedFinder.collisionRegionMax = 150 * u.mm
    # traccc seeding assumes a constant field along z
    seedFinder.bFieldInZ = bFieldInZ
    s.addAlgorithm(
        traccc.TracccDeviceTripletSeedingAlg(
            level=logLevel,
            detector=detector,
            inputSpacepoints="traccc_device_spacepoints",
            outputSeeds="traccc_device_seeds",
            seedFinder=seedFinder,
            seedFilter=traccc.SeedFilterConfig(),
        )
    )

    # example config settings for trk param
    s.addAlgorithm(
        traccc.TracccDeviceTrackParamEstimationAlg(
            level=logLevel,
            detector=detector,
            magneticField=magneticField,
            inputMeasurements="traccc_device_measurements",
            inputSpacepoints="traccc_device_spacepoints",
            inputSeeds="traccc_device_seeds",
            outputTrackParameters="traccc_device_track_parameters",
            estimation=traccc.TrackParamsEstimationConfig(),
        )
    )

    # example config settings for track finding
    finding = traccc.FindingConfig()
    finding.min_track_candidates_per_track = 6
    finding.max_track_candidates_per_track = 100
    finding.max_num_tracks_per_measurement = 1
    finding.min_measurement_voting_fraction = 1.0
    s.addAlgorithm(
        traccc.TracccDeviceTrackFindingAlg(
            level=logLevel,
            detector=detector,
            magneticField=magneticField,
            inputMeasurements="traccc_device_measurements",
            inputTrackParameters="traccc_device_track_parameters",
            outputTracks="traccc_device_tracks",
            finding=finding,
        )
    )

    # Copy the tracks and the (sorted) measurements they refer to back to the
    # host and convert them into ACTS tracks
    s.addAlgorithm(
        traccc.TracccDeviceToHostMeasurementCopyAlg(
            level=logLevel,
            copyTool=copyTool,
            input="traccc_device_measurements",
            output="traccc_host_measurements_sorted",
        )
    )

    s.addAlgorithm(
        traccc.TracccDeviceToHostTrackCopyAlg(
            level=logLevel,
            copyTool=copyTool,
            input="traccc_device_tracks",
            output="traccc_host_tracks",
            inputMeasurements="traccc_host_measurements_sorted",
        )
    )

    s.addAlgorithm(
        traccc.TracccTrackToActsAlg(
            level=logLevel,
            detector=detector,
            inputTracks="traccc_host_tracks",
            inputMeasurements="traccc_host_measurements_sorted",
            outputTracks="traccc_tracks",
        )
    )


def main():
    parser = argparse.ArgumentParser(
        description="ODD simulation + traccc device reconstruction"
    )
    parser.add_argument("--output", "-o", type=Path, default=Path.cwd() / "traccc_odd")
    parser.add_argument("--events", "-n", type=int, default=10)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--backend", choices=["cuda"], default="cuda")
    parser.add_argument("--gun-particles", type=int, default=4)
    parser.add_argument("--gun-multiplicity", type=int, default=50)
    parser.add_argument(
        "--gun-pt-range",
        nargs=2,
        type=float,
        default=[1.0, 10.0],
        metavar=("PT_MIN", "PT_MAX"),
    )
    parser.add_argument(
        "--gun-eta-range",
        nargs=2,
        type=float,
        default=[-3.0, 3.0],
        metavar=("ETA_MIN", "ETA_MAX"),
    )
    parser.add_argument("--ttbar", action="store_true")
    parser.add_argument("--ttbar-pu", type=int, default=200)
    parser.add_argument(
        "--output-root", default=True, action=argparse.BooleanOptionalAction
    )
    parser.add_argument(
        "--material",
        type=Path,
        default=None,
        help="Gen3 (keyed JSON) ODD material map, e.g. from material_mapping.py. "
        "The ROOT map shipped with the ODD is keyed by Gen1 identifiers and "
        "cannot be used. Without it the geometry has no material.",
    )
    parser.add_argument(
        "--log-level",
        choices=["VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
    )
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    logLevel = getattr(acts.logging, args.log_level)

    # ── Geometry ────────────────────────────────────────────────────────────
    geoDir = getOpenDataDetectorDirectory()
    oddMaterialDeco = None
    if args.material is not None:
        oddMaterialDeco = acts.json.MaterialMapDecorator(
            acts.json.TrackingGeometryMaterialJsonConverter().fromFile(
                str(args.material)
            )
        )
    else:
        print("WARNING: no --material given, running without detector material")
    oddDetector = getOpenDataDetector(
        odd_dir=geoDir, materialDecorator=oddMaterialDeco, gen3=True
    )
    trackingGeometry = oddDetector.trackingGeometry()
    gContext = acts.GeometryContext.dangerouslyDefaultConstruct()

    bFieldValue = acts.Vector3(0.0, 0.0, 2.0 * u.T)
    field = acts.ConstantBField(bFieldValue)
    rnd = acts.examples.RandomNumbers(seed=42)

    digiConfig = makeGen3SmearingConfig(
        trackingGeometry, args.output / "odd-digi-smearing-config-gen3.json"
    )

    # ── traccc device setup ─────────────────────────────────────────────────
    detrayDetector, surfaceIds = makeDetrayDetector(
        trackingGeometry, gContext, logLevel
    )

    detectorConfig = traccc.TracccDetector.Config()
    detectorConfig.backend = makeBackend(args.backend, logLevel)
    detectorConfig.trackingGeometry = trackingGeometry
    detectorConfig.detrayDetector = detrayDetector
    detectorConfig.surfaceIds = surfaceIds
    detector = traccc.TracccDetector(detectorConfig, logLevel)

    fieldConfig = traccc.TracccMagneticField.Config()
    fieldConfig.backend = detectorConfig.backend
    fieldConfig.field = traccc.makeTracccHostMagneticField(
        acts.covfie.makeCovfieField(field)
    )
    magneticField = traccc.TracccMagneticField(fieldConfig, logLevel)
    # One copy tool, shared by all host <-> device copy algorithms
    copyTool = traccc.TracccCopyTool(detectorConfig.backend)

    # ── Sequencer ───────────────────────────────────────────────────────────
    # All traccc device algorithms share one stream/queue of the backend
    s = acts.examples.Sequencer(
        events=args.events,
        skip=args.skip,
        numThreads=1,
        logLevel=logLevel,
        outputDir=str(args.output),
        trackFpes=False,
    )

    # ── Simulation ──────────────────────────────────────────────────────────
    if args.ttbar:
        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=args.ttbar_pu,
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
            ),
            rnd=rnd,
        )
        addGenParticleSelection(
            s,
            ParticleSelectorConfig(
                rho=(0.0, 24 * u.mm),
                absZ=(0.0, 1.0 * u.m),
                eta=(-3.0, 3.0),
                pt=(150 * u.MeV, None),
            ),
        )
    else:
        addParticleGun(
            s,
            MomentumConfig(
                args.gun_pt_range[0] * u.GeV,
                args.gun_pt_range[1] * u.GeV,
                transverse=True,
            ),
            EtaConfig(args.gun_eta_range[0], args.gun_eta_range[1]),
            PhiConfig(0.0, 360.0 * u.degree),
            ParticleConfig(
                args.gun_particles, acts.PdgParticle.eMuon, randomizeCharge=True
            ),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns
                ),
            ),
            multiplicity=args.gun_multiplicity,
            rnd=rnd,
        )

    addFatras(s, trackingGeometry, field, enableInteractions=True, rnd=rnd)

    addSimParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            removeNeutral=True,
        ),
    )

    addDigitization(s, trackingGeometry, field, digiConfigFile=digiConfig, rnd=rnd)

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(1.0 * u.GeV, None),
            eta=(-3.0, 3.0),
            measurements=(9, None),
            removeNeutral=True,
        ),
    )

    # ── traccc device reconstruction ────────────────────────────────────────
    addTracccDeviceReconstruction(
        s, detector, magneticField, copyTool, logLevel, bFieldValue[2]
    )

    # ── Truth matching and performance ──────────────────────────────────────
    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=logLevel,
            inputTracks="traccc_tracks",
            inputParticles="particles_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="traccc_track_particle_matching",
            outputParticleTrackMatching="traccc_particle_track_matching",
            doubleMatching=False,
            matchingRatio=0.5,
        )
    )

    if args.output_root:
        from acts.examples.root import (
            RootTrackSummaryWriter,
            RootPatternRecognitionPerformanceWriter,
        )

        s.addWriter(
            RootTrackSummaryWriter(
                level=logLevel,
                inputTracks="traccc_tracks",
                inputParticles="particles_selected",
                inputTrackParticleMatching="traccc_track_particle_matching",
                filePath=str(args.output / "tracksummary_traccc.root"),
                writeCovMat=True,
            )
        )
        s.addWriter(
            RootPatternRecognitionPerformanceWriter(
                level=logLevel,
                inputTracks="traccc_tracks",
                inputParticles="particles_selected",
                inputParticleMeasurementsMap="particle_measurements_map",
                inputTrackParticleMatching="traccc_track_particle_matching",
                inputParticleTrackMatching="traccc_particle_track_matching",
                filePath=str(args.output / "performance_finding_traccc.root"),
            )
        )

    s.run()


if __name__ == "__main__":
    main()
