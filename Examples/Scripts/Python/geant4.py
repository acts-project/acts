#!/usr/bin/env python3
from typing import Optional, Union
from pathlib import Path

import acts
import acts.examples
import acts.examples.dd4hep
from acts.examples.geant4 import Geant4Simulation, geant4SimulationConfig
from acts.examples.geant4.dd4hep import DDG4DetectorConstruction
from common import getOpenDataDetector
from fatras import addFatrasWriters

u = acts.UnitConstants


def addGeant4(
    s: acts.examples.Sequencer,
    geometryService: acts.examples.dd4hep.DD4hepGeometryService,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    seed: Optional[int] = None,
    preselectParticles: bool = True,
) -> acts.examples.Sequencer:
    """This function steers the detector simulation using Geant4

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Geant4 steps (returned from addGeant4)
    trackingGeometry : tracking geometry
    field : magnetic field
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    seed : int, None
        random number generator seed
    """

    if int(s.config.logLevel) <= int(acts.logging.DEBUG):
        acts.examples.dump_args_calls(locals())

    # Selector
    if preselectParticles:
        particles_selected = "particles_selected"
        s.addAlgorithm(
            acts.examples.ParticleSelector(
                level=s.config.logLevel,
                inputParticles="particles_input",
                outputParticles=particles_selected,
            )
        )
    else:
        particles_selected = "particles_input"

    g4detector = DDG4DetectorConstruction(geometryService)
    g4conf = geant4SimulationConfig(
        level=s.config.logLevel,
        detector=g4detector,
        inputParticles="particles_input",
        trackingGeometry=trackingGeometry,
        magneticField=field,
    )
    g4conf.outputSimHits = "simhits"
    g4conf.outputParticlesInitial = "particles_initial"
    g4conf.outputParticlesFinal = "particles_final"
    g4conf.seed = seed

    # Simulation
    alg = Geant4Simulation(
        level=s.config.logLevel,
        config=g4conf,
    )

    # Sequencer
    s.addAlgorithm(alg)

    # Output
    addFatrasWriters(s, None, None)

    return s


def runGeant4(
    geometryService,
    trackingGeometry,
    field,
    outputDir,
    s: acts.examples.Sequencer = None,
):
    from particle_gun import addParticleGun, EtaConfig

    s = s or acts.examples.Sequencer(events=100, numThreads=1)
    s.config.logLevel = acts.logging.INFO
    rnd = acts.examples.RandomNumbers()
    seed = 42
    s = addParticleGun(
        s,
        EtaConfig(-2.0, 2.0),
        rnd=rnd,
    )
    outputDir = Path(outputDir)
    return addGeant4(
        s,
        geometryService,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        seed=seed,
    )


if "__main__" == __name__:
    detector, trackingGeometry, decorators = getOpenDataDetector()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runGeant4(detector.geometryService, trackingGeometry, field, Path.cwd()).run()
