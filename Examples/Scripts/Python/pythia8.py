#!/usr/bin/env python3
from typing import Optional, Union
from pathlib import Path
from collections.abc import Iterable, Sequence

import acts
import acts.examples

u = acts.UnitConstants


def addPythia8(
    s: acts.examples.Sequencer,
    rnd: Optional[acts.examples.RandomNumbers] = None,
    nhard: int = 1,
    npileup: int = 200,
    beam: Optional[
        Union[acts.PdgParticle, Sequence[acts.PdgParticle]]
    ] = None,  # default: acts.PdgParticle.eProton
    cmsEnergy: Optional[float] = None,  # default: 14 * acts.UnitConstants.TeV
    hardProcess: Optional[Sequence[str]] = None,  # default: ["HardQCD:all = on"]
    pileupProcess: Sequence[str] = ["SoftQCD:all = on"],
    vtxGen: Optional[acts.examples.EventGenerator.VertexGenerator] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    printParticles: bool = False,
    returnEvGen: bool = False,
) -> acts.examples.Sequencer:
    """This function steers the particle generation using Pythia8

    NB. this is a reimplementation of common.addPythia8, which is maintained for now for compatibility.

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the particle gun steps (returned from addParticleGun)
    rnd : RandomNumbers, None
        random number generator
    nhard, npileup : int, 1, 200
        Number of hard-scatter and pileup vertices
    beam : PdgParticle|[PdgParticle,PdgParticle], eProton
        beam particle(s)
    cmsEnergy : float, 14 TeV
        CMS energy
    hardProcess, pileupProcess : [str], ["HardQCD:all = on"], ["SoftQCD:all = on"]
        hard and pileup processes
    vtxGen : VertexGenerator, None
        vertex generator module
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    printParticles : bool, False
        print generated particles
    """

    if int(s.config.logLevel) <= int(acts.logging.DEBUG):
        acts.examples.dump_args_calls(locals())

    # Preliminaries
    rnd = rnd or acts.examples.RandomNumbers()
    vtxGen = vtxGen or acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
    )
    if not isinstance(beam, Iterable):
        beam = (beam, beam)

    generators = []
    if nhard is not None and nhard > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=nhard),
                vertex=vtxGen,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=s.config.logLevel,
                    **acts.examples.defaultKWArgs(
                        pdgBeam0=beam[0],
                        pdgBeam1=beam[1],
                        cmsEnergy=cmsEnergy,
                        settings=hardProcess,
                    ),
                ),
            )
        )
    if npileup > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=npileup),
                vertex=vtxGen,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=s.config.logLevel,
                    **acts.examples.defaultKWArgs(
                        pdgBeam0=beam[0],
                        pdgBeam1=beam[1],
                        cmsEnergy=cmsEnergy,
                        settings=pileupProcess,
                    ),
                ),
            )
        )

    # Input
    evGen = acts.examples.EventGenerator(
        level=s.config.logLevel,
        generators=generators,
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    if printParticles:
        s.addAlgorithm(
            acts.examples.ParticlesPrinter(
                level=s.config.logLevel, inputParticles=evGen.config.outputParticles
            )
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        s.addWriter(
            acts.examples.CsvParticleWriter(
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
            acts.examples.RootParticleWriter(
                level=s.config.logLevel,
                inputParticles=evGen.config.outputParticles,
                filePath=str(outputDirRoot / "pythia8_particles.root"),
            )
        )

    return evGen if returnEvGen else s


def runPythia8(
    outputDir,
    outputRoot: bool = True,
    outputCsv: bool = True,
    s: acts.examples.Sequencer = None,
):
    # Preliminaries
    rnd = acts.examples.RandomNumbers()
    outputDir = Path(outputDir)

    # Sequencer
    s = s or acts.examples.Sequencer(
        events=10, numThreads=-1, logLevel=acts.logging.INFO
    )

    return addPythia8(
        s,
        rnd=rnd,
        outputDirCsv=outputDir / "csv" if outputCsv else None,
        outputDirRoot=outputDir if outputRoot else None,
    )


if "__main__" == __name__:
    runPythia8(Path.cwd()).run()
