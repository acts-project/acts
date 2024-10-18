#!/usr/bin/env python3

import os
import argparse
import pathlib
import math

import acts
import acts.examples
import acts.examples.traccc
from acts import GeometryContext
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    ParticleSelectorConfig,
)
from acts.examples.dd4hep import (
    DD4hepDetector,
    DD4hepDetectorOptions,
    DD4hepGeometryService,
)
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

# def create_gen1_tracking_geometry():
#     return getOpenDataDetector()

# def create_gen2_tracking_geometry():

if __name__ == "__main__":
    u = acts.UnitConstants

    parser = argparse.ArgumentParser(description="Full chain with the OpenDataDetector")
    parser.add_argument(
        "--output",
        "-o",
        help="Output directory",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "odd_output",
    )
    parser.add_argument(
        "--events", "-n", help="Number of events", type=int, default=100
    )
    parser.add_argument(
        "--gun-particles",
        help="Multiplicity (no. of particles) of the particle gun",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--gun-multiplicity",
        help="Multiplicity (no. of vertices) of the particle gun",
        type=int,
        default=200,
    )
    parser.add_argument(
        "--gun-eta-range",
        nargs=2,
        help="Eta range of the particle gun",
        type=float,
        default=[-3.0, 3.0],
    )
    parser.add_argument(
        "--gun-pt-range",
        nargs=2,
        help="Pt range of the particle gun (GeV)",
        type=float,
        default=[1.0 * u.GeV, 10.0 * u.GeV],
    )
    parser.add_argument(
        "--digi-config", help="Digitization configuration file", type=pathlib.Path
    )
    parser.add_argument(
        "--material-config", help="Material map configuration file", type=pathlib.Path
    )
    parser.add_argument(
        "--ambi-solver",
        help="Set which ambiguity solver to use, default is the classical one",
        type=str,
        choices=["greedy", "scoring", "ML"],
        default="greedy",
    )
    parser.add_argument(
        "--ambi-config",
        help="Set the configuration file for the Score Based ambiguity resolution",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "ambi_config.json",
    )

    parser.add_argument(
        "--MLSeedFilter",
        help="Use the Ml seed filter to select seed after the seeding step",
        action="store_true",
    )
    parser.add_argument(
        "--reco",
        help="Switch reco on/off",
        default=True,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--output-root",
        help="Switch root output on/off",
        default=True,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--output-csv",
        help="Switch csv output on/off",
        default=True,
        action=argparse.BooleanOptionalAction,
    )

    args = parser.parse_args()

    outputDir = args.output
    geoDir = getOpenDataDetectorDirectory()

    oddMaterialMap = (
        args.material_config
        if args.material_config
        else geoDir / "data/odd-material-maps.root"
    )

    oddDigiConfig = (
        args.digi_config
        if args.digi_config
        else geoDir / "config/odd-digi-geometric-config.json"
    )

    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

    odd_dir = getOpenDataDetectorDirectory()
    odd_xml = odd_dir / "xml" / "OpenDataDetector.xml"

    volumeRadiusCutsMap = {
        28: [850.0],  # LStrip negative z
        30: [850.0],  # LStrip positive z
        23: [400.0, 550.0],  # SStrip negative z
        25: [400.0, 550.0],  # SStrip positive z
        16: [100.0],  # Pixels negative z
        18: [100.0],  # Pixels positive z
    }

    def geoid_hook(geoid, surface):
        gctx = acts.GeometryContext()
        if geoid.volume() in volumeRadiusCutsMap:
            r = math.sqrt(surface.center(gctx)[0] ** 2 + surface.center(gctx)[1] ** 2)

            geoid.setExtra(1)
            for cut in volumeRadiusCutsMap[geoid.volume()]:
                if r > cut:
                    geoid.setExtra(geoid.extra() + 1)

        return geoid

    dd4hepConfig = acts.examples.dd4hep.DD4hepGeometryService.Config(
        xmlFileNames=[str(odd_xml)],
        logLevel=acts.logging.INFO,
        dd4hepLogLevel=acts.logging.INFO,
        geometryIdentifierHook=acts.GeometryIdentifierHook(geoid_hook),
    )
    detector = acts.examples.dd4hep.DD4hepDetector()

    mdecorator = acts.examples.RootMaterialDecorator(
        fileName=str(odd_dir / "data/odd-material-maps.root"),
        level=acts.logging.INFO,
    )

    trackingGeometry, decorators = detector.finalize(dd4hepConfig, mdecorator)

    dd4hepIdGeoIdMap = acts.examples.dd4hep.createDD4hepIdGeoIdMap(trackingGeometry)
    dd4hepIdGeoIdValueMap = {}
    for key, value in dd4hepIdGeoIdMap.items():
        dd4hepIdGeoIdValueMap[key] = value

    geoContext = acts.GeometryContext()
    cOptions = DD4hepDetectorOptions(logLevel=acts.logging.INFO, emulateToGraph="")
    acts.examples.dd4hep.attachDD4hepGeoIdMapper(cOptions, dd4hepIdGeoIdValueMap)
    [final_detector, contextors, store] = detector.finalize(geoContext, cOptions)

    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(
        events=args.events,
        numThreads=1,
        outputDir=str(outputDir),
    )

    addParticleGun(
        s,
        MomentumConfig(
            args.gun_pt_range[0] * u.GeV,
            args.gun_pt_range[1] * u.GeV,
            transverse=True,
        ),
        EtaConfig(args.gun_eta_range[0], args.gun_eta_range[1]),
        PhiConfig(0.0, 360.0 * u.degree),
        ParticleConfig(
            args.gun_particles, acts.PdgParticle.eMuon, randomizeCharge=True
        ),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns),
        ),
        multiplicity=args.gun_multiplicity,
        rnd=rnd,
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        preSelectParticles=(ParticleSelectorConfig()),
        enableInteractions=True,
        outputDirRoot=outputDir if args.output_root else None,
        outputDirCsv=outputDir if args.output_csv else None,
        rnd=rnd,
    )

    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(
            str(oddDigiConfig),
        ),
        surfaceByIdentifier=trackingGeometry.geoIdSurfaceMap(),
        randomNumbers=acts.examples.RandomNumbers(),
        inputSimHits="simhits",
        outputMeasurements="measurements",
        doMerge=True,
        doOutputCells=True,
        doClusterization=True,
    )

    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.VERBOSE)

    s.addAlgorithm(digiAlg)

    detrayStoreOptions = acts.ActsPythonBindings.detray.DetrayConverter.Options()
    detrayStoreOptions.convertMaterial = False
    detrayStore = acts.examples.traccc.convertDetectorHost(
        geoContext, final_detector, detrayStoreOptions
    )

    recoAlg = acts.examples.traccc.ReconstructionChainHostAlgorithm(
        level=acts.logging.INFO,
        field=field,
        detrayStore=detrayStore,
        trackingGeometry=trackingGeometry,
        inputCells="cells",
        inputMeasurements="measurements",
        outputTracks="tracks",
        digitizationConfigs=acts.examples.readDigiConfigFromJson(
            str(oddDigiConfig),
        ),
    )

    s.addAlgorithm(recoAlg)

    s.run()
