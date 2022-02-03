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
    defaults=[1 * u.GeV, 10 * u.GeV, True],
)
EtaConfig = namedtuple(
    "EtaConfig", ["min", "max", "uniform"], defaults=[-4.0, 4.0, True]
)
ParticleConfig = namedtuple(
    "ParticleConfig",
    ["num", "pdg", "randomizeCharge"],
    defaults=[10, PdgParticle.eMuon, True],
)


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

    print("momentumConfig=", momentumConfig)
    print("etaConfig=", etaConfig)
    print("particleConfig=", particleConfig)
    # Preliminaries
    rnd = rnd or RandomNumbers(seed=228)

    # Input
    if vtxGen is None:
        vtxGen = GaussianVertexGenerator()
        vtxGen.stddev = Vector4(0, 0, 0, 0)

    ptclGen = ParametricParticleGenerator(
        p=(momentumConfig.min, momentumConfig.max),
        pTransverse=momentumConfig.transverse,
        eta=(etaConfig.min, etaConfig.max),
        etaUniform=etaConfig.uniform,
        numParticles=particleConfig.num,
        pdg=particleConfig.pdg,
        randomizeCharge=particleConfig.randomizeCharge,
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
        momentumConfig=MomentumConfig(1 * u.GeV, 10 * u.GeV, transverse=False),
        etaConfig=EtaConfig(-4.0, 4.0, uniform=False),
        particleConfig=ParticleConfig(
            2, pdg=acts.PdgParticle.eMuon, randomizeCharge=False
        ),
        printParticles=True,
    )


if "__main__" == __name__:
    import os

    runParticleGun(os.getcwd()).run()
