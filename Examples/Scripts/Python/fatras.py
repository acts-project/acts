#!/usr/bin/env python3
from typing import Optional, Union
from pathlib import Path

import acts
import acts.examples

u = acts.UnitConstants


def addFatras(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    rnd: Optional[acts.examples.RandomNumbers] = None,
) -> acts.examples.Sequencer:
    """This function steers the detector simulation using Fatras

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Fatras steps (returned from addFatras)
    trackingGeometry : tracking geometry
    field : magnetic field
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    rnd : RandomNumbers, None
        random number generator
    """

    # Preliminaries
    rnd = rnd or acts.examples.RandomNumbers()

    # Selector
    selector = acts.examples.ParticleSelector(
        level=s.config.logLevel,
        inputParticles="particles_input",
        outputParticles="particles_selected",
    )

    # Simulation
    alg = acts.examples.FatrasSimulation(
        level=s.config.logLevel,
        inputParticles=selector.config.outputParticles,
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
    )

    # Sequencer
    s.addAlgorithm(selector)
    s.addAlgorithm(alg)

    # Output
    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()
        s.addWriter(
            acts.examples.CsvParticleWriter(
                level=s.config.logLevel,
                outputDir=str(outputDirCsv),
                inputParticles="particles_final",
                outputStem="particles_final",
            )
        )

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        s.addWriter(
            acts.examples.RootParticleWriter(
                level=s.config.logLevel,
                inputParticles="particles_final",
                filePath=str(outputDirRoot / "fatras_particles_final.root"),
            )
        )

    if outputDirCsv is not None:
        s.addWriter(
            acts.examples.CsvParticleWriter(
                level=s.config.logLevel,
                outputDir=str(outputDirCsv),
                inputParticles="particles_initial",
                outputStem="particles_initial",
            )
        )

    if outputDirRoot is not None:
        s.addWriter(
            acts.examples.RootParticleWriter(
                level=s.config.logLevel,
                inputParticles="particles_initial",
                filePath=str(outputDirRoot / "fatras_particles_initial.root"),
            )
        )

    if outputDirCsv is not None:
        s.addWriter(
            acts.examples.CsvSimHitWriter(
                level=s.config.logLevel,
                inputSimHits=alg.config.outputSimHits,
                outputDir=str(outputDirCsv),
                outputStem="hits",
            )
        )

    if outputDirRoot is not None:
        s.addWriter(
            acts.examples.RootSimHitWriter(
                level=s.config.logLevel,
                inputSimHits=alg.config.outputSimHits,
                filePath=str(outputDirRoot / "hits.root"),
            )
        )

    return s


def runFatras(trackingGeometry, field, outputDir, s: acts.examples.Sequencer = None):
    from particle_gun import addParticleGun, EtaConfig

    s = s or acts.examples.Sequencer(events=100, numThreads=-1)
    s.config.logLevel = acts.logging.INFO
    rnd = acts.examples.RandomNumbers()
    s = addParticleGun(
        s,
        EtaConfig(-2.0, 2.0),
        rnd=rnd,
    )
    outputDir = Path(outputDir)
    return addFatras(
        s,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
    )


if "__main__" == __name__:
    import os

    gdc = acts.examples.GenericDetector.Config()
    detector = acts.examples.GenericDetector()
    trackingGeometry, contextDecorators = detector.finalize(gdc, None)

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runFatras(trackingGeometry, field, os.getcwd()).run()
