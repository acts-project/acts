#!/usr/bin/env python3
from pathlib import Path
import matplotlib.pyplot as plt

import acts
import acts.examples
from acts.examples.simulation import addParticleGun, addFatras, addDigitization, EtaConfig, ParticleConfig, PhiConfig, MomentumConfig, addDigiParticleSelection, ParticleSelectorConfig
from acts.examples import DispatchMeasurements, DispatchTrack, DispatchParticles, DispatchParticleMeasurementsMap

u = acts.UnitConstants

def dispatchFunction(dispatchMeasurements: DispatchMeasurements, 
              dispatchParticles: DispatchParticles, 
              dispatchParticleMeasurementsMap: DispatchParticleMeasurementsMap) -> list[DispatchTrack]:
    dispatchTracks = []
        
    print("Received measurements:")
    # make a list of (x,y,z) coordinates of the measurements
    for x, y, z in zip(dispatchMeasurements.x, dispatchMeasurements.y, dispatchMeasurements.z):
        print ((x, y, z))
    
    return dispatchTracks

def runDispatching(trackingGeometry, field, outputDir, digiConfigFile, s: acts.examples.Sequencer = None):
    # Sequencer 
    s = s or acts.examples.Sequencer(events=1000, numThreads=1)
    s.config.logLevel = acts.logging.INFO
    rnd = acts.examples.RandomNumbers()
    addParticleGun(
            s,
            MomentumConfig(1.0 * u.GeV, 2.0 * u.GeV, transverse=True),
            EtaConfig(-0.5, 0.5, uniform=True),
            PhiConfig(0.0, 360.0 * u.degree),
            ParticleConfig(10, acts.PdgParticle.eMuon, randomizeCharge=True),
            multiplicity=1,
            rnd=rnd,
        )  
    outputDir = Path(outputDir)
    
    # Add the dispatching algorithm to the sequencer
    addFatras(
        s,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
    )
        
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )
    
    dpaConfig = acts.examples.PatternDispatchAlgorithm.Config()
    dpaConfig.dispatchFunction = dispatchFunction
    dpaConfig.trackingGeometry = trackingGeometry
    dpaConfig.inputMeasurements = "measurements"
    dpaConfig.outputProtoTracks = "protoTracks"
    
    dpa = acts.examples.PatternDispatchAlgorithm(dpaConfig, acts.logging.INFO)
    s.addAlgorithm(dpa)
    

    return s


if "__main__" == __name__:
    gdc = acts.examples.GenericDetector.Config()
    gd = acts.examples.GenericDetector(gdc)
    trackingGeometry = gd.trackingGeometry()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    digiConfigFile = (
     "generic-digi-smearing-config.json"
    )
    
    runDispatching(trackingGeometry, field, digiConfigFile=digiConfigFile, outputDir=Path.cwd()).run()
