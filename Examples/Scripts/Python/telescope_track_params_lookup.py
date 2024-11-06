#!/usr/bin/env python3

import argparse

import acts
import acts.examples
from acts.examples.simulation import (
    addParticleGun,
    addFatras,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
)

u = acts.UnitConstants

def estimateLookup(
    trackingGeometry,
    numEvents,
    outputPath):

    # Set up the magnetic field
    field = acts.ConstantBField(acts.Vector3(50 * u.T, 0, 0))

    # Fatras simulation of muons
    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(events=numEvents, numThreads=1, logLevel=acts.logging.INFO)

    vertexGen=acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 0, 0),
        mean=acts.Vector4(0, 9, 0, 0)
    )

    addParticleGun(
        s=s,
        etaConfig=EtaConfig(10.0, 10.0),
        phiConfig=PhiConfig(0, 0),
        momentumConfig=MomentumConfig(0.5 * u.GeV, 10 * u.GeV),
        particleConfig=ParticleConfig(1, acts.PdgParticle.eMuon, False),
        multiplicity=1,
        rnd=rnd,
        vtxGen = vertexGen
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        inputParticles="particles_input",
        outputSimHits="sim_hits",
        rnd=rnd,
        preSelectParticles=None
    )

    # Set up the track lookup grid writer
    jsonWriterConfig = acts.examples.JsonTrackParamsLookupWriter.Config(
        path=outputPath
    )
    jsonWriter = acts.examples.JsonTrackParamsLookupWriter(jsonWriterConfig)

    # Set up the track estimation algorithm
    surfaces = list(trackingGeometry.geoIdSurfaceMap().values())
    refSurface = surfaces[0]
    refGeometryId = refSurface.geometryId()

    trackEstConfig = acts.examples.TrackParamsLookupEstimation.Config(
        refLayers={refGeometryId : refSurface},
        bins=(1, 1000),
        inputHits="sim_hits",
        inputParticles="particles_input",
        trackLookupGridWriters = [jsonWriter]
    )
    trackEstAlg = acts.examples.TrackParamsLookupEstimation(trackEstConfig, acts.logging.INFO)

    s.addAlgorithm(trackEstAlg)

    s.run()

def validateLookup(
    trackingGeometry,
    numEvents,
    inputPath):

    # Set up the magnetic field
    field = acts.ConstantBField(acts.Vector3(50 * u.T, 0, 0))

    # Fatras simulation of muons
    rnd = acts.examples.RandomNumbers(seed=92)

    s = acts.examples.Sequencer(events=numEvents, numThreads=1, logLevel=acts.logging.INFO)

    vertexGen=acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 0, 0),
        mean=acts.Vector4(0, 9, 0, 0)
    )

    addParticleGun(
        s=s,
        etaConfig=EtaConfig(10.0, 10.0),
        phiConfig=PhiConfig(0, 0),
        momentumConfig=MomentumConfig(0.5 * u.GeV, 10 * u.GeV),
        particleConfig=ParticleConfig(1, acts.PdgParticle.eMuon, False),
        multiplicity=1,
        rnd=rnd,
        vtxGen = vertexGen
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        inputParticles="particles_input",
        outputSimHits="sim_hits",
        rnd=rnd,
        preSelectParticles=None
    )

    # Set up the track lookup grid reader
    surfaces = list(trackingGeometry.geoIdSurfaceMap().values())
    refSurface = surfaces[0]
    refGeometryId = refSurface.geometryId()

    jsonReaderConfig = acts.examples.JsonTrackParamsLookupReader.Config(
        refLayers={refGeometryId : refSurface},
        bins=(1, 1000))
    jsonReader = acts.examples.JsonTrackParamsLookupReader(jsonReaderConfig)

    lookupConfig = acts.examples.TrackParamsLookupProvider.Config(
        jsonReader,
        inputPath)
    lookup = acts.examples.TrackParamsLookupProvider(lookupConfig)

    refSurf = list(trackingGeometry.geoIdSurfaceMap().values())[0]

    validaterConfig = acts.examples.TrackParamsLookupValidation.Config()
    validaterConfig.refLayers={refSurf.geometryId(): refSurf}
    validaterConfig.lookup=lookup
    validaterConfig.inputHits="sim_hits"
    validaterConfig.inputParticles="particles_input"
    validaterConfig.outputIpPars="ip_pars"
    validaterConfig.outputRefLayerPars="ref_layer_pars"
    validaterConfig.outputIpParsEst="ip_pars_est"
    validaterConfig.outputRefLayerParsEst="ref_layer_pars_est"
    alg = acts.examples.TrackParamsLookupValidation(validaterConfig, acts.logging.INFO)

    s.addAlgorithm(alg)

    # writer
    writerConfig = acts.examples.RootTrackParamsValidationWriter.Config()

    writerConfig.inputIpPars="ip_pars"
    writerConfig.inputRefLayerPars="ref_layer_pars"
    writerConfig.inputIpParsEst="ip_pars_est"
    writerConfig.inputRefLayerParsEst="ref_layer_pars_est"
    writerConfig.path=inputPath.replace(".json", ".root")
    writerConfig.treeName="test"

    writer = acts.examples.RootTrackParamsValidationWriter(writerConfig, acts.logging.INFO)

    s.addWriter(writer)

    s.run()

if __name__ == "__main__":
    p = argparse.ArgumentParser()

    p.add_argument(
        "-nest", "--estimate", type=int, default=10000000, help="Number of events for lookup estimation"
    )
    p.add_argument(
        "-nval", "--validate", type=int, default=100000, help="Number of events for lookup validation"
    )
    p.add_argument(
        "-o", "--output", type=str, default="/home/romanurmanov/tools/acts/acts_telescope_geo/ActsTelescopeGeometryDevelopment_build/test2", help="Output lookup file name"
    )

    args = p.parse_args()

    # Initialize the geometry
    detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
        bounds=[4, 10],
        positions=[30, 60, 90],
        stereos=[0, 0, 0],
        binValue=2,
        surfaceType=0
    )

    estimateLookup(
        trackingGeometry,
        args.estimate,
        args.output)
    
    validateLookup(
        trackingGeometry,
        args.validate,
        args.output + ".json")
    