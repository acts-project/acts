"""Measure Geant4 transport, optionally followed by ODD track reconstruction."""

import argparse
import json
from pathlib import Path
import resource
import sys
from time import perf_counter

import acts
import acts.examples
from acts.examples.simulation import (
    EtaConfig,
    MomentumConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    addDigiParticleSelection,
    addDigitization,
    addGeant4,
    addParticleGun,
    addSimWriters,
)


class RunTiming(acts.examples.IAlgorithm):
    """Sample the boundaries of the Geant4 run from the sequencer thread."""

    def __init__(self):
        super().__init__("RunTiming", acts.logging.WARNING)

    def initialize(self):
        self.initialized = perf_counter()
        return acts.examples.ProcessCode.SUCCESS

    def execute(self, context):
        return acts.examples.ProcessCode.SUCCESS

    def finalize(self):
        self.finalized = perf_counter()
        return acts.examples.ProcessCode.SUCCESS


def add_reconstruction(sequencer, geometry, field, rng):
    from acts.examples.reconstruction import (
        CkfConfig,
        TrackSelectorConfig,
        addCKFTracks,
        addSeeding,
    )

    configs = Path(__file__).resolve().parents[2] / "Configs"
    u = acts.UnitConstants
    addDigitization(
        sequencer,
        geometry,
        field,
        digiConfigFile=configs / "odd-digi-smearing-config.json",
        rnd=rng,
    )
    addDigiParticleSelection(
        sequencer,
        ParticleSelectorConfig(
            pt=(1 * u.GeV, None), measurements=(9, None), removeNeutral=True
        ),
    )
    addSeeding(
        sequencer,
        geometry,
        field,
        geoSelectionConfigFile=configs / "odd-seeding-config.json",
        initialSigmas=[u.mm, u.mm, u.degree, u.degree, 0.0, u.ns],
        initialSigmaQoverPt=0.1 * u.e / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1.0] * 6,
        particleHypothesis=acts.ParticleHypothesis.muon,
    )
    addCKFTracks(
        sequencer,
        geometry,
        field,
        TrackSelectorConfig(nMeasurementsMin=7, maxHoles=2, maxOutliers=2),
        CkfConfig(
            chi2CutOffMeasurement=15.0,
            chi2CutOffOutlier=25.0,
            numMeasurementsCutOff=2,
            seedDeduplication=True,
            stayOnSeed=True,
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--geometry", choices=("telescope", "odd"), required=True)
    parser.add_argument("--field", choices=("constant", "mapped"), required=True)
    parser.add_argument("--events", type=int, default=20)
    parser.add_argument("--particles", type=int, default=1)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--acts-threads", type=int, default=1)
    parser.add_argument(
        "--g4-workers",
        type=int,
        default=0,
        help="Use native Geant4 workers (0 selects serial transport)",
    )
    parser.add_argument(
        "--reconstruct",
        action="store_true",
        help="Run ODD digitization, seeding, CKF and track truth matching",
    )
    args = parser.parse_args()
    if min(args.events, args.particles) < 1 or args.skip < 0:
        parser.error(
            "events and particles must be positive and skip must be nonnegative"
        )
    if args.acts_threads < 1 or args.g4_workers < 0:
        parser.error(
            "ACTS threads must be positive and Geant4 workers must be nonnegative"
        )
    if args.acts_threads != 1 and args.g4_workers == 0:
        parser.error("Concurrent ACTS events require native Geant4 workers")
    if args.reconstruct and args.geometry != "odd":
        parser.error("Track reconstruction requires the ODD geometry")
    args.output.mkdir(parents=True, exist_ok=False)

    setup_start = perf_counter()
    level = acts.logging.WARNING
    u = acts.UnitConstants
    if args.geometry == "telescope":
        detector = acts.examples.TelescopeDetector(
            positions=[30, 60, 90, 120], stereos=[0] * 4, bounds=[200, 200]
        )
    else:
        from acts.examples.odd import getOpenDataDetector

        detector = getOpenDataDetector(logLevel=level, buildCalorimeter=False)

    if args.field == "constant":
        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    else:
        solenoid = acts.SolenoidBField(
            radius=1.2 * u.m, length=6 * u.m, bMagCenter=2 * u.T, nCoils=100
        )
        # Offset the 500 mm grid from ODD sensor boundaries. At a shared
        # boundary, bin lookup can round a step just below it into the next bin.
        z_offset = 0.125 * u.mm
        field = acts.solenoidFieldMap(
            rlim=(0, 50 * u.m),
            zlim=(-50.25 * u.m + z_offset, 50.25 * u.m + z_offset),
            nbins=(101, 202),
            field=solenoid,
        )

    sequencer = acts.examples.Sequencer(
        events=args.events,
        skip=args.skip,
        numThreads=args.acts_threads,
        trackFpes=False,
        outputDir=str(args.output),
        logLevel=level,
    )
    rng = acts.examples.RandomNumbers(seed=args.seed)
    addParticleGun(
        sequencer,
        momentumConfig=MomentumConfig(5 * u.GeV, 5 * u.GeV),
        etaConfig=EtaConfig(1.0, 1.5),
        particleConfig=ParticleConfig(args.particles, acts.PdgParticle.eMuon, False),
        rnd=rng,
    )
    addGeant4(
        sequencer,
        detector,
        detector.trackingGeometry(),
        field,
        rnd=rng,
        physicsList="FTFP_BERT",
        propagatorLargestAcceptableStep=1 * u.m,
        **({"numGeant4Workers": args.g4_workers} if args.g4_workers else {}),
    )
    if args.reconstruct:
        add_reconstruction(sequencer, detector.trackingGeometry(), field, rng)
    # Initialization and finalization run in registration order. This marker
    # brackets event processing and end-run before the ROOT files are closed.
    timing = RunTiming()
    sequencer.addAlgorithm(timing)
    addSimWriters(sequencer, field=field, outputDirRoot=args.output)
    if args.reconstruct:
        from acts.examples.reconstruction import addTrackWriters

        addTrackWriters(
            sequencer, name="ckf", outputDirRoot=args.output, writeCovMat=True
        )
    run_start = perf_counter()
    sequencer.run()
    run_end = perf_counter()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    metrics = {
        "geometry": args.geometry,
        "field": args.field,
        "events": args.events,
        "acts_threads": args.acts_threads,
        "g4_workers": args.g4_workers,
        "reconstruct": args.reconstruct,
        "particles_per_event": args.particles,
        "skip": args.skip,
        "seed": args.seed,
        "setup_s": run_start - setup_start,
        "initialize_s": timing.initialized - run_start,
        "event_window_s": timing.finalized - timing.initialized,
        "run_s": run_end - run_start,
        "peak_rss_mib": rss / (1024**2 if sys.platform == "darwin" else 1024),
    }
    (args.output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")


if __name__ == "__main__":
    main()
