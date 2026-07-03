#!/usr/bin/env python3

import os
import json
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon, Circle
from matplotlib.collections import PatchCollection, LineCollection
from pathlib import Path

import acts
import acts.examples
from acts.json import MaterialMapJsonConverter, TrackingGeometryJsonConverter
from acts.examples.odd import getOpenDataDetector
from acts.examples import (
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
)

from acts.examples.json import (
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
)

def polyArea(face):
    x = [item[0] for item in face]
    y = [item[1] for item in face]
    # shoelace formula for computing area of polygons
    return 0.5 * abs(np.dot(x, np.roll(y,-1)) - np.dot(y, np.roll(x,-1)))
   
def plot2D(vis, projection = "rz", outpath="rzProjection.png"):
    N = len(vis.surfaces)
    poly_patches = [Polygon(face, closed=True) for face in vis.surfaces]     # face = [[x1,y1], [x2,y2],...]   
    poly_collection = PatchCollection(poly_patches, alpha=0.5)
    poly_collection.set_array(vis.color.rgb)
    print(vis.color.rgb)
    #print([item[0] for item in vis.surfaces[0]])
    #print(vis.surfaces[0])
    # Plotting part 
    fig, ax1 = plt.subplots(1,1)
    ax1.add_collection(poly_collection)
    
    if projection == "rz":
        ax1.set_xlabel('z')
        ax1.set_ylabel('r')

    if projection == "xy":
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')

    
    for face in vis.surfaces:
        if polyArea(face) < 1e-8:
            ax1.plot([item[0] for item in face], [item[1] for item in face], color=(vis.color.rgb[0]/255, vis.color.rgb[1]/255, vis.color.rgb[2]/255, 0.5), lw=1)


    plt.savefig(outpath, dpi=300)

def runGeometry(
    trackingGeometry,
    decorators,
    outputDir: Path,
    events=1,
    outputPy=True,
    outputCsv=False,
    outputSurfacesJson=False,
    serializeGeometryJson=False,
):
    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0
        ithread = 0

        context = AlgorithmContext(ialg, ievt, eventStore, ithread)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            if not os.path.isdir(outputDir / "csv"):
                os.makedirs(outputDir / "csv")
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(outputDir / "csv"),
                writePerEvent=True,
            )
            writer.write(context)

        if outputPy:
            vis = acts.PyVisualization(projection='rz')
            trackingGeometry.visualize(
                vis,
                context.geoContext,
                portalViewConfig=acts.ViewConfig(visible=False),
                sensitiveViewConfig=acts.ViewConfig(visible=True, color=acts.Color(255,0,0)),
                viewConfig=acts.ViewConfig(visible=False),
            )
            
            plot2D(vis, projection="rz", outpath="rzProjection.png")
           

        if outputSurfacesJson:
            # if not os.path.isdir(outputDir / "json"):
            #    os.makedirs(outputDir / "json")
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(outputDir / "json"),
                writePerEvent=True,
                writeSensitive=True,
            )
            writer.write(context)

            if outputMaterialMap and ievt == 0:
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
                    fileName=str(outputDir / "geometry-map"),
                    writeFormat=JsonFormat.Json,
                )

                jmw.write(trackingGeometry)

        if serializeGeometryJson:
            converter = TrackingGeometryJsonConverter(level=acts.logging.INFO)
            jsonStr = converter.toJson(context.geoContext, trackingGeometry)
            outPath = outputDir / "json" / "tracking-geometry.json"
            outPath.write_text(jsonStr)


if "__main__" == __name__:
    #detector = acts.examples.GenericDetector()
    detector = getOpenDataDetector(gen3=True)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    #print(Path.cwd())
    runGeometry(trackingGeometry, decorators, outputDir=Path.cwd())

    # Uncomment if you want to create the geometry id mapping for DD4hep
    # dd4hepIdGeoIdMap = acts.examples.dd4hep.createDD4hepIdGeoIdMap(trackingGeometry)
    # dd4hepIdGeoIdValueMap = {}
    # for key, value in dd4hepIdGeoIdMap.items():
    #     dd4hepIdGeoIdValueMap[key] = value.value

    # with open('odd-dd4hep-geoid-mapping.json', 'w') as outfile:
    #    json.dump(dd4hepIdGeoIdValueMap, outfile)
