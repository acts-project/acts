#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

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

    s = s or acts.examples.Sequencer(
        events=10000, numThreads=-1, logLevel=acts.logging.INFO
    )
    rnd = acts.examples.RandomNumbers(seed=42)

    if particlesInput is None:
        addParticleGun(
            s,
            EtaConfig(-3.5, 3.5, uniform=True),
            ParticleConfig(4, acts.PdgParticle.eMuon, True),
            PhiConfig(0.0, 360.0 * u.degree),
            multiplicity=1,
            rnd=rnd,
        )
    else:
        # Read input from input collection (e.g. Pythia8 output)
        evGen = acts.examples.RootParticleReader(
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

    # ODD
    from acts.examples.odd import getOpenDataDetector

    detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()

    #detector = acts.examples.GenericDetector()
    #trackingGeometry = detector.trackingGeometry()

    digiConfigFile = (
       Path.cwd()
        / "odd-geo-digi-config.json"
    )
    assert digiConfigFile.exists()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


    runDigitization(trackingGeometry, field, outputDir=Path.cwd(), digiConfigFile = digiConfigFile).run()