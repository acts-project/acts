#!/usr/bin/env python3
"""
Full chain: ODD simulation → digitization → traccc GPU reconstruction
"""

import pathlib
import argparse
import json

import acts
import acts.examples

from acts.examples.simulation import (
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addParticleGun,
    addFatras,
    addDigitization,
)
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

from acts.examples.traccc import (
    TracccSeqAlgorithm,
    ActsSpToTracccAlg,
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
        description="ODD simulation + traccc GPU reconstruction"
    )
    parser.add_argument(
        "--output",
        "-o",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "traccc_full_output",
    )
    parser.add_argument("--events", "-n", type=int, default=10)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--gun-particles", type=int, default=4)
    parser.add_argument("--gun-multiplicity", type=int, default=50)
    parser.add_argument(
        "--gun-pt-range",
        nargs=2,
        type=float,
        default=[1.0, 10.0],
        metavar=("PT_MIN", "PT_MAX"),
    )
    parser.add_argument(
        "--gun-eta-range",
        nargs=2,
        type=float,
        default=[-3.0, 3.0],
        metavar=("ETA_MIN", "ETA_MAX"),
    )
    parser.add_argument("--detector-file", type=pathlib.Path, required=True)
    parser.add_argument("--material-file", type=pathlib.Path, default=None)
    parser.add_argument("--grid-file", type=pathlib.Path, default=None)
    parser.add_argument("--material-map", type=pathlib.Path, default=None)
    parser.add_argument("--digitization-file", type=pathlib.Path, required=True)
    parser.add_argument("--bfield-file", type=pathlib.Path, required=True)
    parser.add_argument("--conditions-file", type=pathlib.Path, default=True)

    parser.add_argument(
        "--output-root", default=True, action=argparse.BooleanOptionalAction
    )
    parser.add_argument(
        "--output-csv",
        help="Switch csv output on/off",
        default=True,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--log-level",
        choices=["VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
    )
    parser.add_argument("--ttbar", action="store_true")
    parser.add_argument("--ttbar-pu", type=int, default=200)
    parser.add_argument("--geant4", action="store_true")
    parser.add_argument("--do-cpu", action="store_true")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    logLevel = getattr(acts.logging, args.log_level)

    # ── Geometry ──────────────────────────────────────────────────────────────
    geoDir = getOpenDataDetectorDirectory()
    actsDir = pathlib.Path(__file__).parent.parent.parent.parent

    oddMaterialMap = geoDir / "data/odd-material-maps.root"
    oddDigiConfig = actsDir / "Examples/Configs/odd-digi-smearing-config.json"

    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)
    detector = getOpenDataDetector(odd_dir=geoDir, materialDecorator=oddMaterialDeco)
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
    rnd = acts.examples.RandomNumbers(seed=42)

    # ── Sequencer ─────────────────────────────────────────────────────────────
    s = acts.examples.Sequencer(
        events=args.events,
        skip=args.skip,
        numThreads=1,
        logLevel=logLevel,
        outputDir=str(args.output),
        trackFpes=False,
    )

    # ── Step 1: Simulation ────────────────────────────────────────────────────
    if args.ttbar:
        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=args.ttbar_pu,
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
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
            ParticleConfig(
                args.gun_particles, acts.PdgParticle.eMuon, randomizeCharge=True
            ),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns
                ),
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

    # ── Step 2: Digitization to Acts measurements ────────────────
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

    # ── Step 3: Convert Acts measurements to traccc format ────────────────────
    m2tCfg = ActsMeasToTracccAlg.Config()
    m2tCfg.inputActsMeasurements = "measurements"
    m2tCfg.pixelVolumes = [16, 17, 18]  # ODD Gen1 pixel volume IDs
    m2tCfg.detectorFile = str(args.detector_file)
    m2tCfg.outputDetrayToActsMap = "detray-to-acts-map"
    m2tCfg.trackingGeometry = trackingGeometry
    m2tCfg.outputTracccMeasurements = "acts-traccc-measurements"
    s.addAlgorithm(ActsMeasToTracccAlg(m2tCfg, logLevel))

    geoSelectionConfigFile = actsDir / "Examples/Configs/odd-seeding-config.json"
    stripGeoSelectionConfigFile = ""

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
    sp2tCfg.inputSpacePoints = "spacepoints"
    sp2tCfg.outputTracccSpacepoints = "acts-traccc-spacepoints"
    s.addAlgorithm(acts.examples.traccc.ActsSpToTracccAlg(sp2tCfg, logLevel))

    # ── Step 3.5: Write CSV files for Traccc standalone comparison ────────────────
    if args.output_csv:
        s.addWriter(
            acts.examples.CsvSpacePointWriter(
                level=logLevel,
                inputSpacePoints="spacepoints",
                outputDir=str(args.output),
            )
        )
        s.addWriter(
            acts.examples.CsvMeasurementWriter(
                level=logLevel,
                inputMeasurements="measurements",
                inputClusters="clusters",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                outputDir=str(args.output),
            )
        )

    # ── Step 4: Run traccc GPU chain ──────────────────────────────────────────
    seqCfg = TracccSeqAlgorithm.Config()
    seqCfg.detectorFile = str(args.detector_file)
    seqCfg.digitizationFile = str(args.digitization_file)
    seqCfg.conditionsFile = str(args.conditions_file)
    seqCfg.materialFile = str(args.material_file or pathlib.Path())
    seqCfg.gridFile = str(args.grid_file or pathlib.Path())
    seqCfg.bfieldFile = str(args.bfield_file)
    seqCfg.inputMeasurements = "acts-traccc-measurements"
    seqCfg.inputSpacepoints = "acts-traccc-spacepoints"
    seqCfg.backend = (
        TracccSeqAlgorithm.Backend.CPU
        if args.do_cpu
        else TracccSeqAlgorithm.Backend.CUDA
    )
    s.addAlgorithm(TracccSeqAlgorithm(seqCfg, logLevel))

    s.run()


if __name__ == "__main__":
    main()
