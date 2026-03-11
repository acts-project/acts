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
from acts.examples.edm4hep import (
    EDM4hepParticleOutputConverter,
    PodioWriter,
    PodioReader,
    EDM4hepSimInputConverter,
)

u = acts.UnitConstants

from acts.examples.simulation import (
    addParticleGun,
    ParticleConfig,
    EtaConfig,
    PhiConfig,
    MomentumConfig,
    addFatras,
    addPythia8,
    TruthJetConfig,
    addTruthJetAlg,
    addDigitization,
    addSimParticleSelection,
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

from acts.examples.root import (
    RootTrackStatesWriter,
    RootTrackSummaryWriter,
    RootTrackFitterPerformanceWriter,
)


def make_sequencer(
    s: acts.examples.Sequencer,
    inputDir: Path,
    outputDir: Path,
    detector,
    digiConfig,
    trackingGeometry,
    args,
):

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    rnd = acts.examples.RandomNumbers(seed=42)

    if args.edm4hep:
        s.addReader(
            PodioReader(
                level=acts.logging.INFO,
                inputPath=str(inputDir / "test.edm4hep.root"),
                outputFrame="events",
                category="events",
            )
        )

        edm4hepReader = acts.examples.edm4hep.EDM4hepSimInputConverter(
            inputFrame="events",
            inputSimHits=[
                "PixelBarrelReadout",
                "PixelEndcapReadout",
                "ShortStripBarrelReadout",
                "ShortStripEndcapReadout",
                "LongStripBarrelReadout",
                "LongStripEndcapReadout",
            ],
            outputParticlesGenerator="particles_generated",
            outputParticlesSimulation="particles_simulated",
            outputSimHits="simhits",
            outputSimVertices="vertices_truth",
            dd4hepDetector=detector,
            trackingGeometry=trackingGeometry,
            sortSimHitsInTime=False,
            particleRMax=1080 * u.mm,
            particleZ=(-3030 * u.mm, 3030 * u.mm),
            particlePtMin=150 * u.MeV,
            level=acts.logging.INFO,
        )
        s.addAlgorithm(edm4hepReader)

        s.addWhiteboardAlias(
            "particles", edm4hepReader.config.outputParticlesSimulation
        )

        addSimParticleSelection(
            s,
            ParticleSelectorConfig(
                rho=(0.0, 24 * u.mm),
                absZ=(0.0, 1.0 * u.m),
                eta=(-3.0, 3.0),
                removeNeutral=True,
            ),
        )
    else:
        addPythia8(
            s,
            nhard=args.hardscatter,
            npileup=args.pileup,
            # hardProcess=[
            #     "HardQCD:all = off",
            #     "HardQCD:gg2bbbar = on",
            #     "HardQCD:qqbar2bbbar = on",
            #     "PartonLevel:ISR = off",
            #     "PartonLevel:FSR = off",
            #     "PartonLevel:MPI = off",
            #     "HadronLevel:all = on",
            #     "511:mayDecay = off",
            #     "521:mayDecay = off",
            #     "531:mayDecay = off",
            #     "541:mayDecay = off",
            #     "5122:mayDecay = off",
            # ],
            hardProcess=[
                "Top:qqbar2ttbar=on",
                "HadronLevel:Decay=off",
                "StringFlav:probQQtoQ = 0.0",
            ],
            # hardProcess=["WeakSingleBoson:ffbar2gmZ = on","SoftQCD:nonDiffractive=on","HardQCD:hardbbbar=on","HadronLevel:all = on"],
            # hardProcess=["Top:qqbar2ttbar=on", "HadronLevel:Decay = off","StringFlav:probQQtoQ = 0.0","ParticleDecays:limitTau0=off"],
            # hardProcess=["WeakSingleBoson:ffbar2gmZ = on","23:onMode = off", "23:onIfAny = 5","ParticleDecays:limitTau0=off"],
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
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
                pt=(150 * u.MeV, None),
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
            pt=(0.150 * u.GeV, None),
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
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        particleHypothesis=acts.ParticleHypothesis.muon,
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

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.AMVF,
        outputDirRoot=outputDir,
        logLevel=acts.logging.FATAL,
    )

    addTruthJetAlg(
        s,
        TruthJetConfig(
            inputTruthParticles="particles_generated",
            inputTracks="tracks",
            doTrackJetMatching=True,
            outputJets="output_jets",
            jetPtMin=10 * u.GeV,
        ),
        loglevel=acts.logging.INFO,
    )

    s.addWriter(
        acts.examples.root.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputJets="output_jets",
            writeJets=True,
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf.root"),
        )
    )


def make_geometry():
    from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

    geoDir = getOpenDataDetectorDirectory()

    oddMaterialMap = geoDir / "data/odd-material-maps.root"
    oddDigiConfig = actsDir / "Examples/Configs/odd-digi-smearing-config.json"
    oddSeedingSel = actsDir / "Examples/Configs/odd-seeding-config.json"
    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

    detector = getOpenDataDetector(odd_dir=geoDir, materialDecorator=oddMaterialDeco)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    return detector, oddDigiConfig, oddSeedingSel, trackingGeometry, decorators


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

        detector, oddDigiConfig, oddSeedingSel, trackingGeometry, decorators = (
            make_geometry()
        )
        make_sequencer(
            s,
            job_out,
            detector,
            trackingGeometry=trackingGeometry,
            digiConfig=oddDigiConfig,
            args=args,
        )

        s.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", "-n", type=int, default=10)
    parser.add_argument("--pileup", "--pu", type=int, default=0)
    parser.add_argument("--hardscatter", "--hs", type=int, default=1)
    parser.add_argument("--threads", "-t", type=int, default=-1)
    parser.add_argument("--procs", type=int, default=1)
    parser.add_argument("--csv", action="store_true")
    parser.add_argument("--edm4hep", action="store_true")
    args = parser.parse_args()

    inputDir = Path.cwd() / "truth_jet_test_input"
    outputDir = Path.cwd() / "trackToTruth_output"

    runs = [int(f.name[1:]) for f in outputDir.glob("r*")]
    next_run = max(max(runs), 0) + 1 if len(runs) > 0 else 1

    outputDir = outputDir / f"r{next_run:03d}"

    print(outputDir)
    outputDir.mkdir(exist_ok=True, parents=True)

    if args.procs == 1:
        s = acts.examples.Sequencer(
            events=args.events,
            numThreads=args.threads,
            logLevel=acts.logging.INFO,
            outputDir=str(outputDir),
            trackFpes=False,
        )

        detector, oddDigiConfig, oddSeedingSel, trackingGeometry, decorators = (
            make_geometry()
        )
        make_sequencer(
            s,
            inputDir,
            outputDir,
            detector,
            trackingGeometry=trackingGeometry,
            digiConfig=oddDigiConfig,
            args=args,
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
