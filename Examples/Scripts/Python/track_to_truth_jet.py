#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import sys
import time

import acts
import acts.examples

actsDir = Path(__file__).parent.parent.parent.parent

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

from acts.examples.simulation import (
    addParticleGun,
    ParticleConfig,
    EtaConfig,
    PhiConfig,
    MomentumConfig,
    TrackToTruthJetConfig,
    addFatras,
    addPythia8,
    addTrackToTruthJetAlg,
    addDigitization,
    ParticleSelectorConfig,
    addDigiParticleSelection,
    addGenParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    addKalmanTracks,
    addVertexFitting,
    VertexFinder,
)


def make_sequencer(
    s: acts.examples.Sequencer, outputDir: Path, detector, digiConfig, geoSel, args
):

    trackingGeometry = detector.trackingGeometry()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    rnd = acts.examples.RandomNumbers(seed=42)

    addPythia8(
        s,
        nhard=args.hardscatter,
        npileup=args.pileup,
        hardProcess=["Top:qqbar2ttbar=on"],
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
        ),
        rnd=rnd,
        outputDirRoot=None,
        outputDirCsv=None,
        writeHepMC3=None,
    )

    # Effective truth level selection for simulation + track reconstruction
    addGenParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            pt=(500 * u.MeV, None),
        ),
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        enableInteractions=True,
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfig,
        rnd=rnd,
        logLevel=acts.logging.ERROR,
    )

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(0.500 * u.GeV, None),
            measurements=(7, None),
            removeNeutral=True,
            removeSecondaries=False,
        ),
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_generated",
        particleHypothesis=acts.ParticleHypothesis.pion,
        seedingAlgorithm=SeedingAlgorithm.TruthEstimated,
        geoSelectionConfigFile=geoSel,
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0 * u.e / u.GeV,
            1 * u.ns,
        ],
        initialSigmaQoverPt=0.1 * u.e / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1.0] * 6,
    )

    reverseFilteringMomThreshold = 0 * u.GeV

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
        reverseFilteringMomThreshold,
        logLevel=acts.logging.FATAL,
    )

    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected-tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=7,
            ),
        )
    )
    s.addWhiteboardAlias("tracks", "selected-tracks")

    s.addAlgorithm(
        acts.examples.ParticleSelector(
            level=acts.logging.INFO,
            inputParticles="particles_generated",
            outputParticles="jet_input_particles",
        )
    )

    truthJetAlg = acts.examples.TruthJetAlgorithm(
        level=acts.logging.DEBUG,
        inputTruthParticles="jet_input_particles",
        outputJets="truth_jets",
        jetPtMin=10 * u.GeV,
    )

    s.addAlgorithm(truthJetAlg)

    # addTrackToTruthJetAlg(
    #    s,
    #    TrackToTruthJetConfig(
    #        inputTracks="tracks",
    #        inputJets="truth_jets",
    #        outputTrackJets="track_jets",
    #        maxDeltaR=0.4,
    #    ),
    #    loglevel=acts.logging.DEBUG,
    # )


def make_geometry():
    from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

    geoDir = getOpenDataDetectorDirectory()
    # acts.examples.dump_args_calls(locals())  # show python binding calls

    oddMaterialMap = geoDir / "data/odd-material-maps.root"
    assert oddMaterialMap.exists(), f"Material map file {oddMaterialMap} does not exist"

    oddDigiConfig = actsDir / "Examples/Configs/odd-digi-smearing-config.json"
    assert oddDigiConfig.exists(), f"Digi config file {oddDigiConfig} does not exist"

    oddSeedingSel = actsDir / "Examples/Configs/odd-seeding-config.json"
    assert (
        oddSeedingSel.exists()
    ), f"Seeding selection file {oddSeedingSel} does not exist"

    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

    detector = getOpenDataDetector(
        odd_dir=geoDir, materialDecorator=oddMaterialDeco, logLevel=acts.logging.INFO
    )

    return detector, oddDigiConfig, oddSeedingSel


def job(index: int, events: int, skip: int, outputDir: Path, args):
    job_out = outputDir / f"proc_{index:>02d}"
    job_out.mkdir(exist_ok=False)

    with (job_out / "out.log").open("w") as log_file:
        os.dup2(log_file.fileno(), sys.stdout.fileno())
        os.dup2(log_file.fileno(), sys.stderr.fileno())

        s = acts.examples.Sequencer(
            events=events,
            skip=skip,
            numThreads=1,
            logLevel=acts.logging.INFO,
            outputDir=str(job_out),
            trackFpes=False,
        )

        detector, oddDigiConfig, oddSeedingSel = make_geometry()
        make_sequencer(
            s,
            job_out,
            detector,
            digiConfig=oddDigiConfig,
            geoSel=oddSeedingSel,
            args=args,
        )

        s.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", "-n", type=int, default=10)
    parser.add_argument("--skip", "-s", type=int, default=0)
    parser.add_argument("--pileup", "--pu", type=int, default=0)
    parser.add_argument("--hardscatter", "--hs", type=int, default=1)
    parser.add_argument("--threads", "-t", type=int, default=-1)
    parser.add_argument("--procs", type=int, default=1)
    parser.add_argument("--csv", action="store_true")
    args = parser.parse_args()

    outputDir = Path.cwd() / "trackToTruth_output"

    runs = [int(f.name[1:]) for f in outputDir.glob("r*")]
    next_run = max(max(runs), 0) + 1 if len(runs) > 0 else 1

    outputDir = outputDir / f"r{next_run:03d}"

    print(outputDir)
    outputDir.mkdir(exist_ok=True, parents=True)

    if args.procs == 1:
        s = acts.examples.Sequencer(
            events=args.events,
            skip=args.skip,
            numThreads=args.threads,
            logLevel=acts.logging.INFO,
            outputDir=str(outputDir),
            trackFpes=False,
        )

        detector, oddDigiConfig, geoSel = make_geometry()
        make_sequencer(
            s, outputDir, detector, digiConfig=oddDigiConfig, geoSel=geoSel, args=args
        )

        s.run()
    else:

        with ProcessPoolExecutor(max_workers=args.procs) as ex:
            futures = []
            per_proc = args.events // args.procs
            for i in range(args.procs):
                skip = i * per_proc
                nevents = per_proc
                if i == args.procs - 1:
                    nevents = args.events - skip

                futures.append(ex.submit(job, i, nevents, skip, outputDir, args))

            spin = r"/-\|"
            i = 0
            while any([not f.done() for f in futures]):
                time.sleep(0.25)
                i += 1
                i = i % len(spin)

                ndone = len([f for f in futures if f.done()])
                sys.stdout.write(f"\r{spin[i]} {ndone} / {len(futures)} done")
            print()

            for i, f in enumerate(as_completed(futures)):
                try:
                    f.result()
                except Exception as e:
                    print(f"Job failed with exception: {e}")


if __name__ == "__main__":
    main()
