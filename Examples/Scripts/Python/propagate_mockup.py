#!/usr/bin/env python3

"""Propagate muons through the GeoModel muon mockup tracking geometry."""

import argparse
from pathlib import Path

import acts

from muon_mockup import buildMuonMockup
from propagation import runPropagation
from acts.examples.simulation import (
    addParticleGun,
    EtaConfig,
    ParticleConfig,
    MomentumConfig,
)
from acts.examples import (
    AlgorithmContext,
    WhiteBoard,
    ObjTrackingGeometryWriter,
)

from acts.examples.root import RootPropagationSummaryWriter, RootPropagationStepsWriter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        default="",
        help="Optional GeoModel SQLite input file",
    )
    parser.add_argument(
        "--mockupDetector",
        choices=["Muon"],
        default="Muon",
        help="Predefined mockup detector",
    )
    parser.add_argument("--outDir", default="./propagation_mockup")
    parser.add_argument("--nEvents", default=10, type=int)
    parser.add_argument(
        "--verboseSteps",
        action="store_true",
        help="Also write detailed propagation steps",
    )
    parser.add_argument(
        "--obj",
        action="store_true",
        help="Dump the mockup geometry into obj visualization",
    )
    args = parser.parse_args()

    output_dir = Path(args.outDir)
    output_dir.mkdir(parents=True, exist_ok=True)
    u = acts.UnitConstants
    field = acts.ConstantBField(acts.Vector3(0, 0, 0 * u.T))
    sterileLogger = not args.verboseSteps

    _, trackingGeometry, _, factoryCache = buildMuonMockup(
        input=args.input,
        mockupDetector=args.mockupDetector,
        logLevel=acts.logging.INFO,
    )

    seq = acts.examples.Sequencer(
        events=args.nEvents,
        numThreads=1,
        outputDir=output_dir,
        logLevel=acts.logging.INFO,
    )
    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        seq,
        ParticleConfig(num=1, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
        EtaConfig(0.5, 1.0),
        MomentumConfig(1 * u.GeV, 100 * u.GeV, transverse=True),
        rnd=rnd,
    )
    trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
        level=acts.logging.WARNING,
        inputParticles="particles_generated",
        outputTrackParameters="params_particles_generated",
    )
    seq.addAlgorithm(trkParamExtractor)

    nav = acts.Navigator(trackingGeometry=trackingGeometry)

    stepper = acts.EigenStepper(field)

    propagator = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

    propagationAlgorithm = acts.examples.PropagationAlgorithm(
        propagatorImpl=propagator,
        level=acts.logging.INFO,
        sterileLogger=sterileLogger,
        inputTrackParameters="params_particles_generated",
        outputSummaryCollection="propagation_summary",
    )
    seq.addAlgorithm(propagationAlgorithm)

    seq.addWriter(
        RootPropagationSummaryWriter(
            level=acts.logging.INFO,
            inputSummaryCollection="propagation_summary",
            filePath=output_dir / "propagation_summary.root",
        )
    )

    if sterileLogger is False:
        seq.addWriter(
            RootPropagationStepsWriter(
                level=acts.logging.INFO,
                collection="propagation_summary",
                filePath=output_dir / "propagation_steps.root",
            )
        )

    if args.obj:
        wb = WhiteBoard(acts.logging.INFO)
        context = AlgorithmContext(0, 0, wb, 10)
        obj_dir = Path(args.outDir) / "obj"
        obj_dir.mkdir(exist_ok=True)
        writer = ObjTrackingGeometryWriter(
            level=acts.logging.INFO, outputDir=str(obj_dir)
        )

        writer.write(context, trackingGeometry)
    seq.run()


if __name__ == "__main__":
    main()
