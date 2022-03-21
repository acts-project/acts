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
from acts import UnitConstants as u, PdgParticle

# Defaults (given as `None` here) use class defaults defined in
# Examples/Algorithms/Generators/ActsExamples/Generators/ParametricParticleGenerator.hpp
MomentumConfig = namedtuple(
    "MomentumConfig",
    ["min", "max", "transverse"],
    defaults=[None, None, None],
)
EtaConfig = namedtuple(
    "EtaConfig", ["min", "max", "uniform"], defaults=[None, None, None]
)
PhiConfig = namedtuple("PhiConfig", ["min", "max"], defaults=[None, None])
ParticleConfig = namedtuple(
    "ParticleConfig",
    ["num", "pdg", "randomizeCharge"],
    defaults=[None, None, None],
)


@acts.examples.NamedTypeArgs(
    momentumConfig=MomentumConfig,
    etaConfig=EtaConfig,
    phiConfig=PhiConfig,
    particleConfig=ParticleConfig,
)
def addParticleGun(
    s: Sequencer,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    momentumConfig: MomentumConfig = MomentumConfig(),
    etaConfig: Optional[EtaConfig] = EtaConfig(),
    phiConfig: Optional[PhiConfig] = PhiConfig(),
    particleConfig: Optional[ParticleConfig] = ParticleConfig(),
    multiplicity: int = 1,
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

    if int(s.config.logLevel) <= int(acts.logging.DEBUG):
        acts.examples.dump_args_calls(locals())

    # Preliminaries
    rnd = rnd or RandomNumbers(seed=228)

    # Input
    evGen = EventGenerator(
        level=s.config.logLevel,
        generators=[
            EventGenerator.Generator(
                multiplicity=FixedMultiplicityGenerator(n=multiplicity),
                vertex=vtxGen
                or acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    **acts.examples.defaultKWArgs(
                        p=(momentumConfig.min, momentumConfig.max),
                        pTransverse=momentumConfig.transverse,
                        eta=(etaConfig.min, etaConfig.max),
                        phi=(phiConfig.min, phiConfig.max),
                        etaUniform=etaConfig.uniform,
                        numParticles=particleConfig.num,
                        pdg=particleConfig.pdg,
                        randomizeCharge=particleConfig.randomizeCharge,
                    )
                ),
            )
        ],
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
        EtaConfig(-4.0, 4.0),
        ParticleConfig(2),
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        printParticles=True,
    )


if "__main__" == __name__:
    runParticleGun(Path.cwd()).run()
