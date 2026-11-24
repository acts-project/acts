"""Measure one serial Geant4 configuration and save its particle and hit output."""

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--geometry", choices=("telescope", "odd"), required=True)
    parser.add_argument("--field", choices=("constant", "mapped"), required=True)
    parser.add_argument("--events", type=int, default=20)
    parser.add_argument("--particles", type=int, default=1)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if min(args.events, args.particles) < 1 or args.skip < 0:
        parser.error(
            "events and particles must be positive and skip must be nonnegative"
        )
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
        numThreads=1,
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
    )
    # Initialization runs in registration order. At this point Geant4 is ready
    # and no output writer has been initialized. Finalization runs in the same
    # order, so this interval includes end-run but excludes closing output files.
    timing = RunTiming()
    sequencer.addAlgorithm(timing)
    addSimWriters(sequencer, field=field, outputDirRoot=args.output)
    run_start = perf_counter()
    sequencer.run()
    run_end = perf_counter()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    metrics = {
        "geometry": args.geometry,
        "field": args.field,
        "events": args.events,
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
