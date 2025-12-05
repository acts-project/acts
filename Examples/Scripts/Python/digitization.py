#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import argparse

import acts
import acts.examples


u = acts.UnitConstants


def runDigitization(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    outputDir: Path,
    digiConfigFile: Path,
    particlesInput: Optional[Path] = None,
    outputRoot: bool = True,
    outputCsv: bool = True,
    s: Optional[acts.examples.Sequencer] = None,
    doMerge: Optional[bool] = None,
) -> acts.examples.Sequencer:
    from acts.examples.simulation import (
        addParticleGun,
        EtaConfig,
        PhiConfig,
        ParticleConfig,
        addFatras,
        addDigitization,
    )
    from acts.examples.root import RootParticleReader

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )
    rnd = acts.examples.RandomNumbers(seed=42)

    if particlesInput is None:
        addParticleGun(
            s,
            EtaConfig(-4.0, 4.0),
            ParticleConfig(4, acts.PdgParticle.eMuon, True),
            PhiConfig(0.0, 360.0 * u.degree),
            multiplicity=2,
            rnd=rnd,
        )
    else:
        # Read input from input collection (e.g. Pythia8 output)
        evGen = RootParticleReader(
            level=s.config.logLevel,
            filePath=str(particlesInput),
            outputParticles="particles_generated",
        )
        s.addReader(evGen)

        s.addWhiteboardAlias(
            "particles_generated_selected", evGen.config.outputParticles
        )

    outputDir = Path(outputDir)
    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        outputDirCsv=outputDir / "csv" if outputCsv else None,
        outputDirRoot=outputDir if outputRoot else None,
        rnd=rnd,
        doMerge=doMerge,
    )

    return s


if "__main__" == __name__:

    # Parse the command line arguments
    p = argparse.ArgumentParser(description="Digitization")
    p.add_argument(
        "--events",
        "-n",
        type=int,
        help="Number of events",
        default=1000,
    )
    p.add_argument(
        "--type",
        "-t",
        type=str,
        help="Type of digitization",
        default="smearing",
        choices=["smearing", "geometric"],
    )
    p.add_argument(
        "--detector",
        "-d",
        type=str,
        help="Output file",
        default="generic",
        choices=["generic", "odd"],
    )
    args = p.parse_args()

    if args.detector == "generic":
        detector = acts.examples.GenericDetector()
    else:
        from acts.examples.odd import getOpenDataDetector

        detector = getOpenDataDetector()

    trackingGeometry = detector.trackingGeometry()

    digiConfigFile = (
        Path(__file__).resolve().parent.parent.parent.parent
        / "Examples/Configs"
        / f"{args.detector}-digi-{args.type}-config.json"
    )
    assert digiConfigFile.exists()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    s = acts.examples.Sequencer(
        events=args.events, numThreads=-1, logLevel=acts.logging.INFO
    )

    runDigitization(
        trackingGeometry,
        field,
        outputDir=Path.cwd(),
        digiConfigFile=digiConfigFile,
        s=s,
    ).run()
