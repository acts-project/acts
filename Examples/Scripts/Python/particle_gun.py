#!/usr/bin/env python3
import os

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
from acts import Vector4, UnitConstants as u


def runParticleGun(
    outputDirCsv = None,
    outputDirRoot = None, 
    pConfig = [1 * u.GeV, 10 * u.GeV, True],
    etaConfig = [-4., 4., True],
    particleConfig = [10, acts.ActsPythonBindings.PdgParticle.eMuon, True],
    vtxGen = None, 
    s=None
    ):
    """ This function steers the particle generation using the particle gun
    
        Parameters
        ----------
        outputDirCsv : str, path, None
            the output folder for the Csv output, None triggers no output
        outputDirRoot : str, path, None
            the output folder for the Root output, None triggers no output
        pConfig : list
            momentum configuration: minimum momentum, maximum momentum, transverse
        etaConfig : list
            pseudorapidity configuration: eta min, eta max, uniform
        particleConfig: list
            partilce configuration: number of particles, particle type, charge flip
        vtxGen : vertex generator, None
            vertex generator module
        s: Sequencer, None
            the sequencer module, can be provided

    """

    s = s or Sequencer(events=10, numThreads=-1)

    # Preliminaries
    rnd = RandomNumbers(seed=228)

    # Input
    if vtxGen is None :
        vtxGen = GaussianVertexGenerator()
        vtxGen.stddev = Vector4(0, 0, 0, 0)

    ptclGen = ParametricParticleGenerator(
        p=(pConfig[0], pConfig[1]), 
        pTransverse = pConfig[2],
        eta=(etaConfig[0], etaConfig[1]), 
        etaUniform = etaConfig[2],
        numParticles=particleConfig[0],
        pdg=particleConfig[1], 
        randomizeCharge=particleConfig[2]
    )

    g = EventGenerator.Generator()
    g.multiplicity = FixedMultiplicityGenerator(n=1)
    g.vertex = vtxGen
    g.particles = ptclGen

    evGen = EventGenerator(
        level=acts.logging.INFO,
        generators=[g],
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    s.addAlgorithm(
        ParticlesPrinter(
            level=acts.logging.INFO, inputParticles=evGen.config.outputParticles
        )
    )

    if outputDirCsv is not None :
        if not os.path.exists(outputDirCsv):
            os.mkdir(outputDirCsv)

        s.addWriter(
            CsvParticleWriter(
                level=acts.logging.INFO,
                inputParticles=evGen.config.outputParticles,
                outputDir=outputDirCsv,
                outputStem="particles",
            )
        )

    if outputDirRoot is not None :
        if not os.path.exists(outputDirRoot):
            os.mkdir(outputDirRoot)

        root_file = outputDirRoot+"/particles.root"

        s.addWriter(
            RootParticleWriter(
                level=acts.logging.INFO,
                inputParticles=evGen.config.outputParticles,
                filePath=root_file,
            )
        )

    return s


if "__main__" == __name__:

    current_dir = str(os.getcwd())    
    csv_dir = current_dir+'/csv'
    root_dir = current_dir+'/root'
    runParticleGun(csv_dir, root_dir).run()
