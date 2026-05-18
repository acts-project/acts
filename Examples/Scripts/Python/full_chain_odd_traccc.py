#!/usr/bin/env python3
"""
Full chain: ODD simulation → digitization → traccc GPU reco → Acts track validation
"""

import pathlib
import argparse
import json

import acts
import acts.examples

# from acts.examples.root import (RootTrackSummaryWriter, RootTrackFinderPerformanceWriter)
from acts.examples.simulation import (
    MomentumConfig, EtaConfig, PhiConfig, ParticleConfig,
    addParticleGun, addFatras, addDigitization,
)
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

from acts.examples.traccc import (
    TracccSeqAlgorithm, ActsSpToTracccAlg,
    ActsMeasToTracccAlg,
)

from acts.examples.simulation import (
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    addParticleGun,
    addPythia8,
    addGenParticleSelection,
    addFatras,
    addGeant4,
    addSimParticleSelection,
    addDigitization,
    addDigiParticleSelection,
)

u = acts.UnitConstants

def main():
    parser = argparse.ArgumentParser(
        description="ODD simulation + traccc GPU reconstruction")
    parser.add_argument("--output", "-o", type=pathlib.Path,
                        default=pathlib.Path.cwd() / "traccc_full_output")
    parser.add_argument("--events", "-n", type=int, default=10)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--gun-particles", type=int, default=4)
    parser.add_argument("--gun-multiplicity", type=int, default=50)
    parser.add_argument("--gun-pt-range", nargs=2, type=float,
                        default=[1.0, 10.0], metavar=("PT_MIN", "PT_MAX"))
    parser.add_argument("--gun-eta-range", nargs=2, type=float,
                        default=[-3.0, 3.0], metavar=("ETA_MIN", "ETA_MAX"))
    parser.add_argument("--detector-file",  type=pathlib.Path, required=True)
    parser.add_argument("--material-file",  type=pathlib.Path, default=None)
    parser.add_argument("--grid-file",      type=pathlib.Path, default=None)
    parser.add_argument("--material-map", type=pathlib.Path, default=None)
    parser.add_argument("--digitization-file", type=pathlib.Path, required=True)
    parser.add_argument("--bfield-file", type=pathlib.Path, required=True)
    parser.add_argument("--conditions-file", type=pathlib.Path, default=True)
    parser.add_argument("--detray-json-dir", type=pathlib.Path, default=None,
                        help="Directory with pre-generated detray JSON files. "
                             "If not set, they will be generated from the Acts geometry.")
    parser.add_argument("--output-root", default=True,
                        action=argparse.BooleanOptionalAction)
    parser.add_argument(
        "--output-csv",
        help="Switch csv output on/off",
        default=True,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument("--log-level",
                        choices=["VERBOSE","DEBUG","INFO","WARNING","ERROR"],
                        default="INFO")
    parser.add_argument("--ttbar", action="store_true")
    parser.add_argument("--ttbar-pu", type=int, default=200)
    parser.add_argument("--geant4", action="store_true")
    parser.add_argument("--do-cpu", dest="do-cpu", action="store_true")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    logLevel = getattr(acts.logging, args.log_level)

    # ── Geometry ──────────────────────────────────────────────────────────────
    geoDir  = getOpenDataDetectorDirectory()
    actsDir = pathlib.Path(__file__).parent.parent.parent.parent

    oddMaterialMap = geoDir / "data/odd-material-maps.root"
    oddDigiConfig = actsDir / "Examples/Configs/odd-digi-smearing-config.json"

    oddMaterialDeco  = acts.IMaterialDecorator.fromFile(oddMaterialMap)
    detector         = getOpenDataDetector(odd_dir=geoDir,
                                           materialDecorator=oddMaterialDeco)
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
    rnd   = acts.examples.RandomNumbers(seed=42)

    # # ── Detray JSON directory ─────────────────────────────────────────────────
    # detrayJsonDir = args.detray_json_dir or (args.output / "detray_json")
    # detrayJsonDir.mkdir(parents=True, exist_ok=True)

    # ── Sequencer ─────────────────────────────────────────────────────────────
    s = acts.examples.Sequencer(
        events=args.events,
        skip=args.skip,
        numThreads=1,
        logLevel=logLevel,
        outputDir=str(args.output),
        trackFpes=False,
    )

    # ── Step 1: Convert Acts geometry to detray JSON and make detray→Acts map ─────────
    # geoProvCfg = detray_examples.ActsToDetrayDetectorAlg.Config()
    # geoProvCfg.trackingGeometry     = trackingGeometry
    # geoProvCfg.beampipeVolumeName   = "BeamPipe"
    # geoProvCfg.outputDetrayToActsMap = "detray-to-acts-map"
    # geoProvCfg.outputJsonDir        = str(detrayJsonDir)
    # s.addAlgorithm(
    #     detray_examples.ActsToDetrayDetectorAlg(geoProvCfg, logLevel))

    # ── Step 2: Simulation ────────────────────────────────────────────────────
    if args.ttbar:
        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=args.ttbar_pu,
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125*u.mm, 0.0125*u.mm, 55.5*u.mm, 5.0*u.ns),
            ),
            rnd=rnd,
            outputDirRoot=args.output if args.output_root else None,
            outputDirCsv=args.output if args.output_csv else None,
        )
        addGenParticleSelection(
            s,
            ParticleSelectorConfig(
                rho=(0.0, 24 * u.mm),
                absZ=(0.0, 1.0 * u.m),
                eta=(-3.0, 3.0),
                pt=(150 * u.MeV, None),
            ),
        )
    else:
        addParticleGun(
            s,
            MomentumConfig(
                args.gun_pt_range[0] * u.GeV,
                args.gun_pt_range[1] * u.GeV,
                transverse=True,
            ),
            EtaConfig(args.gun_eta_range[0], args.gun_eta_range[1]),
            PhiConfig(0.0, 360.0 * u.degree),
            ParticleConfig(args.gun_particles, acts.PdgParticle.eMuon,
                           randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125*u.mm, 0.0125*u.mm, 55.5*u.mm, 1.0*u.ns),
            ),
            multiplicity=args.gun_multiplicity,
            rnd=rnd,
        )

    if args.geant4:
        if s.config.numThreads != 1:
            raise ValueError("Geant4 does not support multi-threading")
        addGeant4(
            s,
            detector,
            trackingGeometry,
            field,
            outputDirRoot=args.output if args.output_root else None,
            outputDirCsv=args.output if args.output_csv else None,
            rnd=rnd,
            killVolume=trackingGeometry.highestTrackingVolume,
            killAfterTime=25 * u.ns,
        )
    else:
        addFatras(
            s,
            trackingGeometry,
            field,
            enableInteractions=True,
            outputDirRoot=args.output if args.output_root else None,
            outputDirCsv=args.output if args.output_csv else None,
            rnd=rnd,
        )

    addSimParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            removeNeutral=True,
        ),
    )


    # ── Step 3: Digitization to Acts measurements ────────────────
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=oddDigiConfig,
        outputDirRoot=args.output if args.output_root else None,
        rnd=rnd,
    )

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(1.0 * u.GeV, None),
            eta=(-3.0, 3.0),
            measurements=(9, None),
            removeNeutral=True,
        ),
    )

    # ── Step 4: Convert Acts measurements to traccc format ────────────────────
    a2tCfg = ActsMeasToTracccAlg.Config()
    a2tCfg.inputActsMeasurements      = "measurements"
    a2tCfg.detectorFile                = str(args.detector_file)
    a2tCfg.outputDetrayToActsMap = "detray-to-acts-map"
    a2tCfg.trackingGeometry = trackingGeometry
    a2tCfg.outputTracccMeasurements   = "acts-traccc-measurements"
    s.addAlgorithm(ActsMeasToTracccAlg(a2tCfg, logLevel))

    geoSelectionConfigFile = actsDir / "Examples/Configs/odd-seeding-config.json"
    stripGeoSelectionConfigFile = ""
    # stripGeoSelectionConfigFile = actsDir / "Examples/Configs/odd-seeding-config.json"

    spAlg = acts.examples.SpacePointMaker(
        level=logLevel,
        inputMeasurements=f"measurements",
        outputSpacePoints=f"spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.json.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
        stripGeometrySelection=(
            acts.examples.json.readJsonGeometryList(str(stripGeoSelectionConfigFile))
            if stripGeoSelectionConfigFile
            else []
        ),
    )
    s.addAlgorithm(spAlg)


    sp2tCfg = acts.examples.traccc.ActsSpToTracccAlg.Config()
    sp2tCfg.inputSpacePoints          = "spacepoints"
    # sp2tCfg.inputActsMeasurements     = "measurements"
    sp2tCfg.outputTracccSpacepoints   = "acts-traccc-spacepoints"
    s.addAlgorithm(acts.examples.traccc.ActsSpToTracccAlg(sp2tCfg, logLevel))


    # # ── Step 4.5: Write CSV files for Traccc standalone comparison ────────────────
    # if args.output_csv:
    #     s.addWriter(
    #         acts.examples.CsvSpacePointWriter(
    #             level=logLevel,
    #             inputSpacePoints="spacepoints",
    #             outputDir=str(args.output),
    #         )
    #     )
    #     s.addWriter(
    #         acts.examples.CsvMeasurementWriter(
    #             level=logLevel,
    #             inputMeasurements="measurements",
    #             inputClusters="clusters",
    #             inputMeasurementSimHitsMap="measurement_simhits_map",
    #             outputDir=str(args.output),
    #         )
    #     )


    # ── Step 5: Run traccc GPU chain ──────────────────────────────────────────
    seqCfg = TracccSeqAlgorithm.Config()
    seqCfg.detectorFile              = str(args.detector_file)
    seqCfg.digitizationFile          = str(args.digitization_file)
    seqCfg.conditionsFile            = str(args.conditions_file)
    seqCfg.materialFile              = str(args.material_file   or pathlib.Path())
    seqCfg.gridFile                  = str(args.grid_file       or pathlib.Path())
    seqCfg.bfieldFile                = str(args.bfield_file)
    seqCfg.inputMeasurements     = "acts-traccc-measurements"
    seqCfg.inputSpacepoints      = "acts-traccc-spacepoints"
    seqCfg.backend = TracccSeqAlgorithm.Backend.CUDA
    # seqCfg.outputTracks          = "traccc-tracks"
    # seqCfg.outputDetrayToActsMap = "detray-to-acts-map-traccc"
    s.addAlgorithm(TracccSeqAlgorithm(seqCfg, logLevel))

    # # ── Step 6: Convert traccc tracks to Acts tracks ───────────────────────────
    # trkCfg = TracccTrackToActsAlg.Config()
    # trkCfg.inputTracccTracks    = "traccc-tracks"
    # trkCfg.inputDetrayToActsMap = "detray-to-acts-map"
    # trkCfg.trackingGeometry     = trackingGeometry
    # trkCfg.outputActsTracks     = "traccc-acts-tracks"
    # s.addAlgorithm(TracccTrackToActsAlg(trkCfg, logLevel))

    # # ── Step 7: Truth matching ────────────────────────────────────────────────
    # matcherCfg = acts.examples.TrackTruthMatcher.Config()
    # matcherCfg.inputTracks                  = "traccc-acts-tracks"
    # matcherCfg.inputParticles               = "particles_simulated"
    # matcherCfg.inputMeasurementParticlesMap = "measurement_particles_map"
    # matcherCfg.outputTrackParticleMatching  = "traccc_track_particle_matching"
    # matcherCfg.outputParticleTrackMatching  = "traccc_particle_track_matching"
    # matcherCfg.doubleMatching               = False
    # matcherCfg.matchingRatio                = 0.5
    # s.addAlgorithm(acts.examples.TrackTruthMatcher(matcherCfg, logLevel))

    # # ── Step 8: Validation ────────────────────────────────────────────
    # if args.output_root:
    #     s.addWriter(
    #         RootTrackSummaryWriter(
    #             level=logLevel,
    #             inputTracks="traccc-acts-tracks",
    #             inputParticles="particles_selected",
    #             inputTrackParticleMatching="traccc_track_particle_matching",
    #             filePath=str(args.output / "traccc_track_summary.root"),
    #             writeCovMat=True,
    #         )
    #     )

    #     s.addWriter(
    #         RootTrackFinderPerformanceWriter(
    #             level=logLevel,
    #             inputTracks="traccc-acts-tracks",
    #             inputParticles="particles_selected",
    #             inputParticleMeasurementsMap="particle_measurements_map",
    #             inputTrackParticleMatching="traccc_track_particle_matching",
    #             inputParticleTrackMatching="traccc_particle_track_matching",
    #             filePath=str(args.output / "traccc_track_finder_performance.root"),
    #         )
    #     )

    # ── Step 9: Do CPU chain + validation ────────────────────────────────────────────
    # if args.do-cpu:
    #     addSeeding(
    #         s,
    #         trackingGeometry,
    #         field,
    #         initialSigmas=[
    #             1 * u.mm,
    #             1 * u.mm,
    #             1 * u.degree,
    #             1 * u.degree,
    #             0 * u.e / u.GeV,
    #             1 * u.ns,
    #         ],
    #         initialSigmaQoverPt=0.1 * u.e / u.GeV,
    #         initialSigmaPtRel=0.1,
    #         initialVarInflation=[1.0] * 6,
    #         particleHypothesis=acts.ParticleHypothesis.muon,
    #         geoSelectionConfigFile=oddSeedingSel,
    #         outputDirRoot=outputDir if args.output_root else None,
    #         outputDirCsv=outputDir if args.output_csv else None,
    #     )
    #     addCKFTracks(
    #         s,
    #         trackingGeometry,
    #         field,
    #         TrackSelectorConfig(
    #             pt=(1.0 * u.GeV if args.ttbar else 0.0, None),
    #             absEta=(None, 3.0),
    #             loc0=(-4.0 * u.mm, 4.0 * u.mm),
    #             nMeasurementsMin=7,
    #             maxHoles=2,
    #             maxOutliers=2,
    #         ),
    #         CkfConfig(
    #             chi2CutOffMeasurement=15.0,
    #             chi2CutOffOutlier=25.0,
    #             numMeasurementsCutOff=2,
    #             seedDeduplication=True,
    #             stayOnSeed=True,
    #             pixelVolumes=[16, 17, 18],
    #             stripVolumes=[23, 24, 25],
    #             maxPixelHoles=1,
    #             maxStripHoles=2,
    #             constrainToVolumes=[
    #                 2,  # beam pipe
    #                 32,
    #                 4,  # beam pip gap
    #                 16,
    #                 17,
    #                 18,  # pixel
    #                 20,  # PST
    #                 23,
    #                 24,
    #                 25,  # short strip
    #                 26,
    #                 8,  # long strip gap
    #                 28,
    #                 29,
    #                 30,  # long strip
    #             ],
    #         ),
    #         outputDirRoot=outputDir if args.output_root else None,
    #         outputDirCsv=outputDir if args.output_csv else None,
    #         writeCovMat=True,
    #     )

    #     addAmbiguityResolution(
    #         s,
    #         AmbiguityResolutionConfig(
    #             maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
    #         ),
    #         outputDirRoot=outputDir if args.output_root else None,
    #         outputDirCsv=outputDir if args.output_csv else None,
    #         writeCovMat=True,
    #     )

    #     if args.output_root:
    #         s.addWriter(
    #             RootTrackSummaryWriter(
    #                 level=logLevel,
    #                 inputTracks="traccc-acts-tracks",
    #                 inputParticles="particles_selected",
    #                 inputTrackParticleMatching="traccc_track_particle_matching",
    #                 filePath=str(args.output / "traccc_track_summary.root"),
    #                 writeCovMat=True,
    #             )
    #         )

    #         s.addWriter(
    #             RootTrackFinderPerformanceWriter(
    #                 level=logLevel,
    #                 inputTracks="traccc-acts-tracks",
    #                 inputParticles="particles_selected",
    #                 inputParticleMeasurementsMap="particle_measurements_map",
    #                 inputTrackParticleMatching="traccc_track_particle_matching",
    #                 inputParticleTrackMatching="traccc_particle_track_matching",
    #                 filePath=str(args.output / "traccc_track_finder_performance.root"),
    #             )
    #         )


    s.run()


if __name__ == "__main__":
    main()