#!/usr/bin/env python3

import argparse
from pathlib import Path

import acts
import acts.examples
import acts.examples.edm4hep
from acts.examples.edm4hep import (
    PodioReader,
    EDM4hepSimInputConverter,
    EDM4hepCaloHitInputConverter,
    CaloCollectionDetectorCodes,
)

u = acts.UnitConstants

ODD_SIM_HIT_COLLECTIONS = [
    "PixelBarrelReadout",
    "PixelEndcapReadout",
    "ShortStripBarrelReadout",
    "ShortStripEndcapReadout",
    "LongStripBarrelReadout",
    "LongStripEndcapReadout",
]

ODD_CALO_HIT_COLLECTIONS = [
    "ECalBarrelCollection",
    "ECalEndcapCollection",
    "HCalBarrelCollection",
    "HCalEndcapCollection",
]

# Calo detector codes for ODD-style EDM4hep collection names (see ArrowCalo defaults).
ODD_CALO_DETECTOR_ENCODER_TABLE = {
    "ECalBarrelCollection": CaloCollectionDetectorCodes.barrel(10),
    "ECalEndcapCollection": CaloCollectionDetectorCodes.endcap(9, 11),
    "HCalBarrelCollection": CaloCollectionDetectorCodes.barrel(13),
    "HCalEndcapCollection": CaloCollectionDetectorCodes.endcap(12, 14),
}


def main():
    parser = argparse.ArgumentParser(
        description="Convert EDM4hep simulation output to parquet tables."
    )
    parser.add_argument(
        "--input", "-i", type=Path, required=True, help="EDM4hep input ROOT file"
    )
    parser.add_argument(
        "--output", "-o", type=Path, required=True, help="Output directory for parquet"
    )
    parser.add_argument("--events", "-n", type=int, default=-1)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--jobs", "-j", type=int, default=-1)
    parser.add_argument(
        "--sim-hits",
        nargs="+",
        default=ODD_SIM_HIT_COLLECTIONS,
        help="EDM4hep SimTrackerHit collection names (defaults to ODD)",
    )
    parser.add_argument(
        "--no-helix",
        action="store_true",
        help="Skip writing perigee d0/z0 (no propagator/field needed)",
    )
    parser.add_argument(
        "--calo-hits",
        nargs="+",
        default=ODD_CALO_HIT_COLLECTIONS,
        help="EDM4hep SimCalorimeterHit collection names (defaults to ODD)",
    )
    parser.add_argument(
        "--calo",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip writing calorimeter cells",
    )
    parser.add_argument(
        "--force", "-f", action="store_true", help="Allow non-empty output directory"
    )
    args = parser.parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input file {args.input} does not exist")

    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()) and not args.force:
        raise FileExistsError(
            f"Output directory {args.output} is not empty; pass --force to proceed"
        )

    try:
        from acts.arrow import particleSchema, caloHitSchema
        from acts.examples.arrow import (
            ArrowParticleOutputConverter,
            ArrowCaloHitOutputConverter,
            ParquetWriter,
        )
    except ImportError as e:
        raise RuntimeError(
            "acts.examples.arrow is not available; "
            "rebuild with ACTS_BUILD_EXAMPLES_PARQUET=ON."
        ) from e

    from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

    geoDir = getOpenDataDetectorDirectory()
    materialMap = geoDir / "data/odd-material-maps.root"
    detector = getOpenDataDetector(
        odd_dir=geoDir, materialDecorator=acts.IMaterialDecorator.fromFile(materialMap)
    )
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))

    s = acts.examples.Sequencer(
        events=args.events if args.events > 0 else None,
        skip=args.skip,
        numThreads=args.jobs,
        outputDir=str(args.output),
        logLevel=acts.logging.INFO,
    )

    s.addReader(
        PodioReader(
            level=acts.logging.INFO,
            inputPath=str(args.input),
            outputFrame="events",
            category="events",
        )
    )

    s.addAlgorithm(
        EDM4hepSimInputConverter(
            level=acts.logging.INFO,
            inputFrame="events",
            inputSimHits=args.sim_hits,
            outputParticlesGenerator="particles_generated",
            outputParticlesSimulation="particles_simulated",
            outputSimHits="simhits",
            outputSimVertices="vertices_truth",
            outputMCParticleMap=("mcparticle_index_map" if args.calo else None),
            dd4hepDetector=detector,
            trackingGeometry=trackingGeometry,
            sortSimHitsInTime=False,
        )
    )

    s.addAlgorithm(
        ArrowParticleOutputConverter(
            level=acts.logging.INFO,
            inputParticles="particles_simulated",
            outputTable="particles_arrow",
            writeHelixParameters=not args.no_helix,
            bField=None if args.no_helix else field,
        )
    )

    parquet_collections = {"particles_arrow": "particles"}
    parquet_schemas = {"particles_arrow": particleSchema()}

    if args.calo:
        s.addAlgorithm(
            EDM4hepCaloHitInputConverter(
                level=acts.logging.INFO,
                inputFrame="events",
                inputCaloHitCollections=args.calo_hits,
                inputMCParticleMap="mcparticle_index_map",
                outputCaloHits="calohits",
                caloDetectorCodesByCollectionName=ODD_CALO_DETECTOR_ENCODER_TABLE,
            )
        )

        s.addAlgorithm(
            ArrowCaloHitOutputConverter(
                level=acts.logging.INFO,
                inputCaloHits="calohits",
                outputTable="calohits_arrow",
            )
        )

        parquet_collections["calohits_arrow"] = "calohits"
        parquet_schemas["calohits_arrow"] = caloHitSchema()

    s.addWriter(
        ParquetWriter(
            level=acts.logging.INFO,
            outputDir=str(args.output),
            collections=parquet_collections,
            expectedSchemas=parquet_schemas,
        )
    )

    s.run()


if __name__ == "__main__":
    main()
