#!/usr/bin/env python3
import os
from typing import Optional, Union, Tuple
from pathlib import Path

from acts.examples import (
    Sequencer,
    RandomNumbers,
    GaussianVertexGenerator,
    EventGenerator,
    FixedMultiplicityGenerator,
    ParametricParticleGenerator,
    CsvParticleWriter,
    ParticlesPrinter,
    RootParticleWriter,
)

import acts
from acts import Vector4, UnitConstants as u, PdgParticle


def addParticleGun(
    s: Sequencer,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    momentumConfig: Tuple[float, float, bool] = (1 * u.GeV, 10 * u.GeV, True),
    etaConfig: Tuple[float, float, bool] = (-4.0, 4.0, True),
    particleConfig: Tuple[int, PdgParticle, bool] = (10, PdgParticle.eMuon, True),
    vtxGen: Optional[EventGenerator.VertexGenerator] = None,
    printParticles: bool = False,
    rnd: Optional[RandomNumbers] = None,
) -> Sequencer:
    """This function steers the particle generation using the particle gun

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the particle gun steps (returned from addParticleGun)
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    momentumConfig : tuple|list
        momentum configuration: minimum momentum, maximum momentum, transverse
    etaConfig : tuple|list
        pseudorapidity configuration: eta min, eta max, uniform
    particleConfig: tuple|list
        partilce configuration: number of particles, particle type, charge flip
    vtxGen : VertexGenerator, None
        vertex generator module
    printParticles : bool, False
        print generated particles
    rnd : RandomNumbers, None
        random number generator
    """

    # Preliminaries
    rnd = rnd or RandomNumbers(seed=228)

    # Input
    if vtxGen is None:
        vtxGen = GaussianVertexGenerator()
        vtxGen.stddev = Vector4(0, 0, 0, 0)

    ptclGen = ParametricParticleGenerator(
        p=(momentumConfig[0], momentumConfig[1]),
        pTransverse=momentumConfig[2],
        eta=(etaConfig[0], etaConfig[1]),
        etaUniform=etaConfig[2],
        numParticles=particleConfig[0],
        pdg=particleConfig[1],
        randomizeCharge=particleConfig[2],
    )

    g = EventGenerator.Generator()
    g.multiplicity = FixedMultiplicityGenerator(n=1)
    g.vertex = vtxGen
    g.particles = ptclGen

    evGen = EventGenerator(
        level=s.config.logLevel,
        generators=[g],
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    if printParticles:
        s.addAlgorithm(
            ParticlesPrinter(
                level=s.config.logLevel, inputParticles=evGen.config.outputParticles
            )
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not os.path.exists(outputDirCsv):
            os.mkdir(outputDirCsv)

        s.addWriter(
            CsvParticleWriter(
                level=s.config.logLevel,
                inputParticles=evGen.config.outputParticles,
                outputDir=str(outputDirCsv),
                outputStem="particles",
            )
        )

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not os.path.exists(outputDirRoot):
            os.mkdir(outputDirRoot)

        s.addWriter(
            RootParticleWriter(
                level=s.config.logLevel,
                inputParticles=evGen.config.outputParticles,
                filePath=str(outputDirRoot / "particles.root"),
            )
        )

    return s


def runParticleGun(outputDir, s=None):
    s = s or Sequencer(events=10, numThreads=-1)
    s.config.logLevel = acts.logging.INFO
    outputDir = Path(outputDir)
    return addParticleGun(
        s,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        momentumConfig=(1 * u.GeV, 10 * u.GeV, False),
        etaConfig=(-4.0, 4.0, False),
        particleConfig=(2, acts.PdgParticle.eMuon, False),
        printParticles=True,
    )


if "__main__" == __name__:
    runParticleGun(os.getcwd()).run()
