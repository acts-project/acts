#!/usr/bin/env python3
from typing import Optional, Union, Tuple
from pathlib import Path
from collections import namedtuple

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

MomentumConfig = namedtuple(
    "MomentumConfig",
    ["min", "max", "transverse"],
    defaults=[None, None, None],
)
EtaConfig = namedtuple(
    "EtaConfig", ["min", "max", "uniform"], defaults=[None, None, None]
)
ParticleConfig = namedtuple(
    "ParticleConfig",
    ["num", "pdg", "randomizeCharge"],
    defaults=[None, None, None],
)


def DefaultKWArgs(**kwargs) -> dict:
    """Removes keyword arguments that are None or a list of all None (eg. [None,None]).
    This keeps the called function's defaults."""
    from collections.abc import Iterable

    return {
        k: v
        for k, v in kwargs.items()
        if not (
            v is None or (isinstance(v, Iterable) and all([vv is None for vv in v]))
        )
    }


def addParticleGun(
    s: Sequencer,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    momentumConfig: MomentumConfig = MomentumConfig(),
    etaConfig: EtaConfig = EtaConfig(),
    particleConfig: ParticleConfig = ParticleConfig(),
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
    momentumConfig : MomentumConfig(min, max, transverse)
        momentum configuration: minimum momentum, maximum momentum, transverse
    etaConfig : EtaConfig(min, max, uniform)
        pseudorapidity configuration: eta min, eta max, uniform
    particleConfig: ParticleConfig(num, pdg, randomizeCharge)
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
        **DefaultKWArgs(
            p=(momentumConfig.min, momentumConfig.max),
            pTransverse=momentumConfig.transverse,
            eta=(etaConfig.min, etaConfig.max),
            etaUniform=etaConfig.uniform,
            numParticles=particleConfig.num,
            pdg=particleConfig.pdg,
            randomizeCharge=particleConfig.randomizeCharge,
        )
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
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

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
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

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
        etaConfig=EtaConfig(-4.0, 4.0),
        particleConfig=ParticleConfig(2),
        printParticles=True,
    )


if "__main__" == __name__:
    import os

    runParticleGun(os.getcwd()).run()
