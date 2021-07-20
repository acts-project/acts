#!/usr/bin/env python3
import os

import acts
import acts.examples

u = acts.UnitConstants


def runPythia8(outputDir, s: acts.examples.Sequencer = None):
    # Preliminaries
    rnd = acts.examples.RandomNumbers()

    nhard = 1
    npileup = 200
    beam0 = acts.PdgParticle.eProton
    beam1 = acts.PdgParticle.eProton
    cmsEnergy = 14 * u.TeV

    vertexGenerator = acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
    )

    generators = []
    if nhard > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=nhard),
                vertex=vertexGenerator,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=acts.logging.INFO,
                    pdgBeam0=beam0,
                    pdgBeam1=beam1,
                    cmsEnergy=cmsEnergy,
                    settings=["HardQCD:all = on"],
                ),
            )
        )
    if npileup > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=npileup),
                vertex=vertexGenerator,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=acts.logging.INFO,
                    pdgBeam0=beam0,
                    pdgBeam1=beam1,
                    cmsEnergy=cmsEnergy,
                    settings=["SoftQCD:all = on"],
                ),
            )
        )

    # Input
    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=generators,
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    # Sequencer
    s = s or acts.examples.Sequencer(
        events=10, numThreads=-1, logLevel=acts.logging.INFO
    )

    s.addReader(evGen)
    s.addWriter(
        acts.examples.RootParticleWriter(
            level=acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            filePath=outputDir + "/pythia8_particles.root",
        )
    )
    s.addWriter(
        acts.examples.CsvParticleWriter(
            level=acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            outputDir=outputDir + "/csv/",
            outputStem="particles",
        )
    )

    return s


if "__main__" == __name__:
    runPythia8(os.getcwd()).run()
