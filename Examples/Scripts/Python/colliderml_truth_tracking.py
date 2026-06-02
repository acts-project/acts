#!/usr/bin/env python3
"""Minimal truth-tracking demonstrator using ColliderML data.

Reads particles and tracker hits from a ColliderML parquet dataset,
converts them to ACTS EDM (SimParticleContainer, SimHitContainer,
MeasurementContainer), selects particles with pT > 1 GeV and >= 5 hits,
runs truth-smeared seeding, and fits tracks with the Kalman filter.

Usage
-----
./run_in_env.sh python3 Examples/Scripts/Python/colliderml_truth_tracking.py \\
    --input  colliderml-sample/CERN__ColliderML-Release-1 \\
    --digi-config build/_deps/odd-src/config/odd-digi-smearing-config.json \\
    --output colliderml_output
"""

import argparse
import pathlib

import acts
import acts.examples
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory
from acts.examples.simulation import (
    addDigiParticleSelection,
    ParticleSelectorConfig,
)
from acts.examples.reconstruction import (
    SeedingAlgorithm,
    TrackSmearingSigmas,
    addSeeding,
    addKalmanTracks,
)

u = acts.UnitConstants


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        "-i",
        type=pathlib.Path,
        required=True,
        help="ColliderML sample root directory "
        "(contains ttbar_pu200_{particles,tracker_hits}/)",
    )
    parser.add_argument(
        "--geo-map",
        type=pathlib.Path,
        default=None,
        help="ColliderML → ACTS geometry ID map CSV "
        "(from generate_colliderml_geo_map.py). "
        "Omit to use direct volume/layer/surface passthrough.",
    )
    parser.add_argument(
        "--digi-config",
        type=pathlib.Path,
        required=True,
        help="ACTS smearing digitisation config JSON "
        "(e.g. odd-digi-smearing-config.json)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "colliderml_output",
        help="Output directory (default: colliderml_output)",
    )
    parser.add_argument(
        "--events",
        "-n",
        type=int,
        default=10,
        help="Number of events to process (default: 10)",
    )
    parser.add_argument(
        "--odd-dir",
        type=pathlib.Path,
        default=None,
        help="ODD XML directory (default: auto-detect)",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel threads (default: 1)",
    )
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Detector + field
    # ------------------------------------------------------------------
    odd_dir = args.odd_dir or getOpenDataDetectorDirectory()
    odd = getOpenDataDetector(odd_dir=odd_dir, withExtraBits=False)
    tgeo = odd.trackingGeometry()
    decorators = odd.contextDecorators()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    # ------------------------------------------------------------------
    # Sequencer
    # ------------------------------------------------------------------
    s = acts.examples.Sequencer(
        events=args.events,
        numThreads=args.jobs,
        logLevel=acts.logging.INFO,
        outputDir=str(args.output),
        # boost::histogram log-axis fills produce benign FPEs in ROOT writers
        failOnUnmaskedFpe=False,
    )

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    # ------------------------------------------------------------------
    # Imports
    # ------------------------------------------------------------------
    from acts.arrow import particleSchema, simHitSchema
    from acts.examples.arrow import (
        loadColliderMLGeoIdMap,
        ColliderMLInputConverter,
        ParquetReader,
    )
    from acts.examples.json import readDigiConfigFromJson
    from acts.examples.root import (
        RootTrackStatesWriter,
        RootTrackSummaryWriter,
        RootTrackFitterPerformanceWriter,
    )

    geo_id_map = loadColliderMLGeoIdMap(args.geo_map) if args.geo_map else {}
    digi_config = readDigiConfigFromJson(str(args.digi_config))

    # ------------------------------------------------------------------
    # ParquetReader: place both tables on the whiteboard per event
    # ------------------------------------------------------------------
    particles_dir = (
        args.input / "ttbar_pu200_particles" / "data" / "ttbar_pu200_particles"
    ).resolve()
    hits_dir = (
        args.input / "ttbar_pu200_tracker_hits" / "data" / "ttbar_pu200_tracker_hits"
    ).resolve()

    s.addReader(
        ParquetReader(
            level=acts.logging.INFO,
            inputDir=str(particles_dir.parent),  # unused but required
            collections={
                "cml_particles": str(particles_dir),
                "cml_hits": str(hits_dir),
            },
            expectedSchemas={
                "cml_particles": particleSchema(),
                "cml_hits": simHitSchema(),
            },
        )
    )

    # ------------------------------------------------------------------
    # ColliderMLInputConverter: particles + simhits + measurements
    # ------------------------------------------------------------------
    s.addAlgorithm(
        ColliderMLInputConverter(
            level=acts.logging.INFO,
            inputParticlesTable="cml_particles",
            inputHitsTable="cml_hits",
            outputParticles="particles",
            outputSimHits="simhits",
            outputMeasurements="measurements",
            outputMeasurementSubset="measurement_subset",
            outputMeasSimHitsMap="measurement_simhits_map",
            outputMeasParticlesMap="measurement_particles_map",
            outputParticleMeasurementsMap="particle_measurements_map",
            trackingGeometry=tgeo,
            digiConfig=digi_config,
            geoIdMap=geo_id_map,
        )
    )

    # ------------------------------------------------------------------
    # Particle selection: pT > 1 GeV, >= 5 measurements, charged only
    # addDigiParticleSelection reads "particles_simulated_selected"
    # ------------------------------------------------------------------
    s.addWhiteboardAlias("particles_simulated_selected", "particles")

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(1.0 * u.GeV, None),
            measurements=(5, None),
            removeNeutral=True,
        ),
    )
    # addDigiParticleSelection aliases output to "particles_selected"

    # ------------------------------------------------------------------
    # Truth-smeared seeding
    # Zero smearing so the KF fits from true parameters
    # ------------------------------------------------------------------
    addSeeding(
        s,
        trackingGeometry=tgeo,
        field=field,
        rnd=rnd,
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        selectedParticles="particles_selected",
        trackSmearingSigmas=TrackSmearingSigmas(
            loc0=0,
            loc0PtA=0,
            loc0PtB=0,
            loc1=0,
            loc1PtA=0,
            loc1PtB=0,
            time=0,
            phi=0,
            theta=0,
            ptRel=0,
        ),
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0 / u.GeV,
            1 * u.ns,
        ],
        initialSigmaQoverPt=0.1 / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1e0, 1e0, 1e0, 1e0, 1e0, 1e0],
        logLevel=acts.logging.INFO,
    )

    # ------------------------------------------------------------------
    # Kalman filter track fitting
    # ------------------------------------------------------------------
    addKalmanTracks(
        s,
        trackingGeometry=tgeo,
        field=field,
        logLevel=acts.logging.INFO,
    )

    # ------------------------------------------------------------------
    # Track selection
    # ------------------------------------------------------------------
    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected_tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=5,
            ),
        )
    )
    s.addWhiteboardAlias("tracks", "selected_tracks")

    # ------------------------------------------------------------------
    # ROOT output
    # ------------------------------------------------------------------
    s.addWriter(
        RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(args.output / "trackstates_kf.root"),
        )
    )
    s.addWriter(
        RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(args.output / "tracksummary_kf.root"),
        )
    )
    s.addWriter(
        RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(args.output / "performance_kf.root"),
        )
    )

    s.run()


if __name__ == "__main__":
    main()
