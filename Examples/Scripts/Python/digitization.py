#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union

import acts
import acts.examples
from acts.examples.odd import getOpenDataDetector
from common import getOpenDataDetectorDirectory

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
        MomentumConfig,
        ParticleConfig,
        addFatras,
        addDigitization,
    )

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )
    rnd = acts.examples.RandomNumbers(seed=42)

    if particlesInput is None:
        addParticleGun(
            s,
            MomentumConfig(0.25 * u.GeV, 100.0 * u.GeV, transverse=True),
            EtaConfig(-3.2, 3.2),
            ParticleConfig(1000, acts.PdgParticle.eMuon, True),
            PhiConfig(0.0, 360.0 * u.degree),
            multiplicity=2,
            rnd=rnd,
        )
    else:
        # Read input from input collection (e.g. Pythia8 output)
        evGen = acts.examples.RootParticleReader(
            level=s.config.logLevel,
            particleCollection="particles_input",
            filePath=str(particlesInput),
            orderedEvents=False,
        )
        s.addReader(evGen)

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

    geoDir = getOpenDataDetectorDirectory()

    detector, trackingGeometry, decorators = getOpenDataDetector(
        geoDir
    )

    #digiConfigFile = geoDir / "config/odd-digi-smearing-config.json"
    #digiConfigFile = geoDir / "config/odd-digi-geometric-config.json"
    #digiConfigFile = geoDir / "config/odd-digi-geometric-new-config.json"
 
    digiConfigFile = Path.cwd() / "odd-digi-geometric-config.json"

    assert digiConfigFile.exists()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runDigitization(trackingGeometry, field, digiConfigFile=digiConfigFile, outputDir=Path.cwd()).run()
