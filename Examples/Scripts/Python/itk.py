#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse

from acts.examples import (
    TGeoDetector,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
    Interval,
)

import acts

from acts import MaterialMapJsonConverter, UnitConstants as u


def runITk(
    trackingGeometry,
    decorators,
    outputDir: Path,
    events=1,
    outputObj=True,
    outputCsv=False,
    outputJson=False,
):

    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0

        context = AlgorithmContext(ialg, ievt, eventStore)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            csv_dir = outputDir / "csv"
            csv_dir.mkdir(exist_ok=True)
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(csv_dir),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            obj_dir = outputDir / "obj"
            obj_dir.mkdir(exist_ok=True)
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO,
                outputDir=str(obj_dir),
            )
            writer.write(context, trackingGeometry)

        if outputJson:
            json_dir = outputDir / "json"
            json_dir.mkdir(exist_ok=True)
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(json_dir),
                writePerEvent=True,
                writeSensitive=True,
            )
            writer.write(context)

            jmConverterCfg = MaterialMapJsonConverter.Config(
                processSensitives=True,
                processApproaches=True,
                processRepresenting=True,
                processBoundaries=True,
                processVolumes=True,
                processNonMaterial=True,
                context=context.geoContext,
            )

            jmw = JsonMaterialWriter(
                level=acts.logging.VERBOSE,
                converterCfg=jmConverterCfg,
                fileName=str(json_dir / "material-map"),
                writeFormat=JsonFormat.Json,
            )

            jmw.write(trackingGeometry)


def buildITkGeometry(geo_dir: Path, material: bool = True):
    Volume = TGeoDetector.Config.Volume
    LayerTriplet = TGeoDetector.Config.LayerTriplet
    equidistant = TGeoDetector.Config.BinningType.equidistant
    arbitrary = TGeoDetector.Config.BinningType.arbitrary

    logger = acts.logging.getLogger("buildITkGeometry")

    matDeco = None
    if material:
        file = geo_dir / "atlas/itk-hgtd/material-maps-ITk-HGTD.json"
        logger.info("Adding material from %s", file.absolute())
        matDeco = acts.IMaterialDecorator.fromFile(
            file,
            level=acts.logging.INFO,
        )

    return TGeoDetector.create(
        fileName=str(geo_dir / "atlas/itk-hgtd/ATLAS-ITk-HGTD.tgeo.root"),
        mdecorator=matDeco,
        buildBeamPipe=True,
        unitScalor=1.0,  # explicit units
        beamPipeRadius=29.0 * u.mm,
        beamPipeHalflengthZ=3000.0 * u.mm,
        beamPipeLayerThickness=0.8 * u.mm,
        surfaceLogLevel=acts.logging.WARNING,
        layerLogLevel=acts.logging.WARNING,
        volumeLogLevel=acts.logging.WARNING,
        volumes=[
            Volume(
                name="InnerPixels",
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet("Pixel::Pixel"),
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                sensitiveNames=LayerTriplet(["Pixel::siLog"]),
                sensitiveAxes=LayerTriplet("YZX"),
                rRange=LayerTriplet((0 * u.mm, 135 * u.mm)),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -250 * u.mm),
                    central=(-250 * u.mm, 250 * u.mm),
                    positive=(250 * u.mm, 3000 * u.mm),
                ),
                splitTolR=LayerTriplet(negative=-1.0, central=5 * u.mm, positive=-1.0),
                splitTolZ=LayerTriplet(
                    negative=10 * u.mm, central=-1.0, positive=10 * u.mm
                ),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(6, equidistant), (10, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(12, equidistant), (6, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
            ),
            Volume(
                name="OuterPixels",
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet("Pixel::Pixel"),
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                sensitiveNames=LayerTriplet(["Pixel::siLog"]),
                sensitiveAxes=LayerTriplet("YZX"),
                rRange=LayerTriplet((135 * u.mm, 350 * u.mm)),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -377 * u.mm),
                    central=(-377 * u.mm, 377 * u.mm),
                    positive=(377 * u.mm, 3000 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    negative=15 * u.mm, central=5 * u.mm, positive=15 * u.mm
                ),
                splitTolZ=LayerTriplet(
                    negative=20 * u.mm, central=-1.0, positive=20 * u.mm
                ),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
            ),
            Volume(
                name="Strips",
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet("*"),
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                sensitiveNames=LayerTriplet(
                    negative=["SCT::ECSensor*"],
                    central=["SCT::BRLSensor*"],
                    positive=["SCT::ECSensor*"],
                ),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    negative=(-1.0, 1050 * u.mm),
                    central=(380 * u.mm, 1050 * u.mm),
                    positive=(-1.0, 1050 * u.mm),
                ),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -1400 * u.mm),
                    central=(-1400 * u.mm, 1400 * u.mm),
                    positive=(1400 * u.mm, 3000 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    negative=-1.0,
                    central=35 * u.mm,
                    positive=-1.0,
                ),
                splitTolZ=LayerTriplet(
                    negative=35 * u.mm, central=-1.0, positive=35 * u.mm
                ),
                binning0=LayerTriplet(
                    negative=[(-1, arbitrary)],
                    central=[(0, equidistant)],
                    positive=[(-1, arbitrary)],
                ),
                binning1=LayerTriplet(
                    negative=[(-1, arbitrary)],
                    central=[(28, equidistant)] * 4,
                    positive=[(-1, arbitrary)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=True,
                barrelMap={"MS": 2, "SS": 4},
                discMap={
                    "EC0": [
                        [384.5, 403.481],
                        [403.481, 427.462],
                        [427.462, 456.442],
                        [456.442, 488.423],
                    ],
                    "EC1": [
                        [489.823, 507.916],
                        [507.916, 535.009],
                        [535.009, 559.101],
                        [559.101, 574.194],
                    ],
                    "EC2": [[575.594, 606.402], [606.402, 637.209]],
                    "EC3": [
                        [638.609, 670.832],
                        [670.832, 697.055],
                        [697.055, 723.278],
                        [723.278, 755.501],
                    ],
                    "EC4": [[756.901, 811.482], [811.482, 866.062]],
                    "EC5": [[867.462, 907.623], [907.623, 967.785]],
                },
            ),
            Volume(
                name="HGTD",
                layers=LayerTriplet(positive=True, central=False, negative=True),
                subVolumeName=LayerTriplet("HGTD::HGTD"),
                binToleranceR=(15 * u.mm, 15 * u.mm),
                binToleranceZ=(15 * u.mm, 15 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                sensitiveNames=LayerTriplet(["HGTD::HGTDSiSensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    negative=(0 * u.mm, 1050 * u.mm),
                    positive=(0 * u.mm, 1050 * u.mm),
                ),
                zRange=LayerTriplet(
                    negative=(-4000 * u.mm, -3000 * u.mm),
                    postive=(3000 * u.mm, 4000 * u.mm),
                ),
                splitTolR=LayerTriplet(-1.0),
                splitTolZ=LayerTriplet(negative=10 * u.mm, positive=10 * u.mm),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
            ),
        ],
    )


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Example script to construct the ITk geometry and write it out to CSV and OBJ formats"
    )
    p.add_argument(
        "geo_dir",
        help="Input directory containing the ITk standalone geometry. Get in touch if you don't have this.",
    )
    p.add_argument(
        "--output-dir",
        default=Path.cwd(),
        type=Path,
        help="Directory to write outputs to",
    )
    p.add_argument(
        "--output-csv", action="store_true", help="Write geometry in CSV format."
    )
    p.add_argument(
        "--output-obj", action="store_true", help="Write geometry in OBJ format."
    )
    p.add_argument(
        "--output-json",
        action="store_true",
        help="Write geometry and material in JSON format.",
    )
    p.add_argument(
        "--no-material", action="store_true", help="Decorate material to the geometry"
    )

    args = p.parse_args()
    args.output_dir.mkdir(exist_ok=True, parents=True)

    geo_example_dir = Path(args.geo_dir)
    assert geo_example_dir.exists(), "Detector example input directory missing"

    detector, trackingGeometry, decorators = buildITkGeometry(
        geo_example_dir,
        material=not args.no_material,
    )

    runITk(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        outputDir=args.output_dir,
        outputCsv=args.output_csv,
        outputObj=args.output_obj,
        outputJson=args.output_json,
    )
