#!/usr/bin/env python3
# This file is part of the ACTS project.
# Copyright (C) 2016 CERN for the benefit of the ACTS project
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Construct DD4hep-backed Gen3, audit surfaces and run straight navigation."""

import argparse
import json
from pathlib import Path

import acts
import acts.examples
from acts.examples.nodd import getNoddDetector
from acts.examples.root import RootPropagationSummaryWriter, RootPropagationStepsWriter
from acts.examples.simulation import (
    addParticleGun,
    EtaConfig,
    MomentumConfig,
    ParticleConfig,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodd-dir", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tracks", type=int, default=128)
    args = parser.parse_args()
    if args.tracks <= 0:
        parser.error("tracks must be positive")
    args.output.mkdir(parents=True, exist_ok=True)
    detector = getNoddDetector(
        nodd_dir=args.nodd_dir, build_dir=args.build_dir, logLevel=acts.logging.WARNING
    )
    report = json.loads(detector.geometryReport())
    (args.output / "surfaces.json").write_text(json.dumps(report))
    sequence = acts.examples.Sequencer(
        events=1, numThreads=1, logLevel=acts.logging.WARNING
    )
    addParticleGun(
        sequence,
        ParticleConfig(num=args.tracks, pdg=acts.PdgParticle.eMuon),
        EtaConfig(-4, 4),
        MomentumConfig(10 * acts.UnitConstants.GeV, 10 * acts.UnitConstants.GeV),
        rnd=acts.examples.RandomNumbers(seed=42),
    )
    sequence.addAlgorithm(
        acts.examples.ParticleTrackParamExtractor(
            level=acts.logging.WARNING,
            inputParticles="particles_generated",
            outputTrackParameters="start_parameters",
        )
    )
    navigator = acts.Navigator(trackingGeometry=detector.trackingGeometry())
    propagator = acts.examples.ConcretePropagator(
        acts.Propagator(acts.StraightLineStepper(), navigator)
    )
    sequence.addAlgorithm(
        acts.examples.PropagationAlgorithm(
            propagatorImpl=propagator,
            level=acts.logging.WARNING,
            sterileLogger=False,
            energyLoss=False,
            multipleScattering=False,
            recordMaterialInteractions=True,
            inputTrackParameters="start_parameters",
            outputSummaryCollection="propagation_summary",
        )
    )
    sequence.addWriter(
        RootPropagationSummaryWriter(
            level=acts.logging.WARNING,
            inputSummaryCollection="propagation_summary",
            filePath=str(args.output / "summary.root"),
        )
    )
    sequence.addWriter(
        RootPropagationStepsWriter(
            level=acts.logging.WARNING,
            collection="propagation_summary",
            filePath=str(args.output / "steps.root"),
        )
    )
    sequence.run()
    print(
        f"Gen3: {len(report['sensors'])} sensors, {len(report['volumes'])} volumes; {args.tracks} straight tracks, seed42"
    )


if __name__ == "__main__":
    main()
