import random
import re

import acts
import argparse
from acts import (
    logging,
    GeometryContext,
    CylindricalContainerBuilder,
    DetectorBuilder,
    GeometryIdGenerator,
)
from acts import geomodel as gm
from acts import examples

from propagation import runPropagation

def merge(a, b):
    if isinstance(a, list) and isinstance(b, list):
        return a + b

    assert isinstance(a, Hierarchy) and isinstance(b, Hierarchy)
    assert a.level == b.level
    assert a.max_level == b.max_level
    assert a.settings == b.settings

    merged = Hierarchy(a.level, a.max_level, a.settings)

    keys = set([*a.keys(), *b.keys()])

    for key in keys:
        if key in a and not key in b:
            merged[key] = a[key]
        elif key in b and not key in a:
            merged[key] = b[key]
        else:
            merged[key] = merge(a[key], b[key])

    return merged


class HierarchySetting:
    '''
    This class will be passed to all child hierarchies by reference (compared to a single bool, that will be passed by value)
    This way, we can manipulate all child trees with one call
    '''
    allow_insertion = True

class Hierarchy(dict):
    '''
    Small utility class, that allows to easily construct fixed depth hierarchies
    After maxlevel is reached, a list is inserted instead of another Hierarchy
    '''
    def __init__(self, level, max_level, settings = HierarchySetting()):
        self.level= level
        self.max_level = max_level
        self.settings = settings

    def __getitem__(self, key):
        if key is None:
            assert self.level < self.max_level
            broadcast = Hierarchy(self.level+1, self.max_level, self.settings)
            for k in self.keys():
                broadcast = merge(broadcast, self[k])

            return broadcast

        if isinstance(key, list):
            assert self.level < self.max_level
            broadcast = Hierarchy(self.level+1, self.max_level, self.settings)
            for k in key:
                broadcast = merge(broadcast, self[k])

            return broadcast

        if self.settings.allow_insertion and not key in self:
            if self.level < self.max_level:
                self[key] = Hierarchy(self.level+1, self.max_level, self.settings)
            else:
                self[key] = []

        return super().__getitem__(key)

    def flatten(self):
        '''
        Flatten out the hierarchy at a given level
        '''
        l = []
        if self.level < self.max_level:
            for key in self.keys():
                l += self[key].flatten()
        else:
            for key in self.keys():
                l += self[key]

        return l



class ItkBuilder:
    def __init__(self, gmFactoryCache, gctx, logLevel):
        self.pyLogger = acts.logging.getLogger("ItkBuilder")
        self.pyLogger.setLevel(logLevel)

        gmSurfaces = [ss[1] for ss in gmFactoryCache.sensitiveSurfaces]
        gmDetElements = [ss[0] for ss in gmFactoryCache.sensitiveSurfaces]

        pixelPattern = r"^barrel_endcap_(-?\d+)_eta_module_(-?\d+)_layer_wheel_(-?\d+)_phi_module_(-?\d+)_side_0_RD53.*$"
        stripPattern = r"^barrel_endcap_(-?\d+)_eta_module_(-?\d+)_layer_wheel_(-?\d+)_phi_module_(-?\d+)_side_(0|1)_split.*$"

        self.index_hierarchy = Hierarchy(0, 4)

        for surfaceIdx in range(len(gmDetElements)):
            detEl = gmDetElements[surfaceIdx]
            if match := re.match(pixelPattern, detEl.databaseEntryName()):
                barrel_endcap = int(match[1])
                hardware = "PIXEL"
                eta = int(match[2])
                layer_wheel =int(match[3])
                phi_module = int(match[4])
            elif match := re.match(stripPattern, detEl.databaseEntryName()):
                barrel_endcap = int(match[1])
                hardware = "STRIP"
                eta = int(match[2])
                layer_wheel =int(match[3])
                phi_module = int(match[4])
            else:
                self.error(f"Could not match {detEl.databaseEntryName()}")
                continue

            self.index_hierarchy[hardware][barrel_endcap][layer_wheel][eta][phi_module].append(gmSurfaces[surfaceIdx])

        self.index_hierarchy.settings.allow_insertion = False
        self.gctx = gctx

        layerCreatorCfg = acts.LayerCreator.Config()
        layerCreatorCfg.surfaceArrayCreator = acts.SurfaceArrayCreator()
        self.layerCreator = acts.LayerCreator(layerCreatorCfg)

        cylVolHelperCfg = acts.CylinderVolumeHelper.Config()
        cylVolHelperCfg.layerArrayCreator = acts.LayerArrayCreator()
        cylVolHelperCfg.trackingVolumeArrayCreator = acts.TrackingVolumeArrayCreator()
        self.cylVolHelper = acts.CylinderVolumeHelper(cylVolHelperCfg, acts.logging.INFO)

        self.rPixInner = 30
        self.rPixOuter = 130
        self.rStripInner = 350
        self.rMax = 1050
        self.zMaxEc = 3000

    def info(self, msg):
        self.pyLogger.log(acts.logging.INFO, msg)

    def error(self, msg):
        self.pyLogger.log(acts.logging.ERROR, msg)

    def buildInnerPixel(self):
        rmin = self.rPixInner
        rmax = self.rPixOuter
        zmax_barrel = 250
        zmax_ec = self.zMaxEc

        self.info("Build inner pixel barrel")
        pixelBarrel = self.index_hierarchy["PIXEL"][0]

        pixel1Vols = []

        pixelBarrelBounds1 = acts.CylinderVolumeBounds(rmin=rmin, rmax=rmax, halfz=zmax_barrel)
        pixelBarrelLayers1 = [ self.layerCreator.cylinderLayer(self.gctx, pixelBarrel[lid].flatten(), 1, 1) for lid in [0,1] ]
        pixel1Vols.append(self.cylVolHelper.createTrackingVolume(self.gctx, pixelBarrelLayers1, pixelBarrelBounds1, 0.0, "PixelBarrelInner"))

        for ec in [-2, 2]:
            self.info(f"Build inner pixel endcaps {ec}")
            pixelEndcap = self.index_hierarchy["PIXEL"][ec]

            sign = ec / abs(ec)
            halfz_ec = 0.5 * (zmax_ec - zmax_barrel)
            shift = sign * (zmax_barrel + halfz_ec)

            pixelEndcapBounds1 = acts.CylinderVolumeBounds(rmin=rmin, rmax=rmax, halfz=halfz_ec)
            pixelEndcapLayers1 = \
                [ self.layerCreator.discLayer(self.gctx, pixelEndcap[[0,2]][eta].flatten(), 1, 1) for eta in range(23) ] + \
                [ self.layerCreator.discLayer(self.gctx, pixelEndcap[1][eta].flatten(), 1, 1) for eta in range(6)]
            pixel1Vols.append(self.cylVolHelper.createTrackingVolume(self.gctx, pixelEndcapLayers1, pixelEndcapBounds1, shift, f"PixelEndcap{ec}Inner"))

        return self.cylVolHelper.createContainerTrackingVolume(self.gctx, pixel1Vols)

    def buildOuterPixel(self):
        pixel2Vols = []

        rmin = self.rPixOuter
        rmax = self.rStripInner
        zmax_barrel = 380
        zmax_ec = self.zMaxEc

        self.info("Build outer pixel barrel")

        pixelBarrelBounds2 = acts.CylinderVolumeBounds(rmin=rmin, rmax=rmax, halfz=zmax_barrel)
        pixelBarrelLayers2 = [ self.layerCreator.cylinderLayer(self.gctx, self.index_hierarchy["PIXEL"][0][lid].flatten(), 1, 1) for lid in [2,3,4] ]
        pixel2Vols.append(self.cylVolHelper.createTrackingVolume(self.gctx, pixelBarrelLayers2, pixelBarrelBounds2, 0.0, "PixelBarrelOuter"))

        for ec in [-2, 2]:
            self.info(f"Build outer pixel endcaps {ec}")
            endcap = self.index_hierarchy["PIXEL"][ec]

            sign = ec / abs(ec)

            rmid1 = 200
            rmid2 = 260

            # Tilted part of endcaps
            tiltedVols = []

            zmax_tilted = 1060
            halfz_tilted = 0.5 * (zmax_tilted - zmax_barrel)
            zshift_tilted = sign * (zmax_barrel + halfz_tilted)

            tiltedLayers1 = [ self.layerCreator.discLayer(self.gctx, endcap[3][eta].flatten(), 1, 1) for eta in range(6) ]
            tiltedBounds1 = acts.CylinderVolumeBounds(rmin=rmin, rmax=rmid1, halfz=halfz_tilted)
            tiltedVols.append(self.cylVolHelper.createTrackingVolume(self.gctx, tiltedLayers1, tiltedBounds1, zshift_tilted,  f"PixelEndcap{ec}Tilted1"))

            tiltedLayers2 = [ self.layerCreator.discLayer(self.gctx, endcap[5][eta].flatten(), 1, 1) for eta in range(8) ]
            tiltedBounds2 = acts.CylinderVolumeBounds(rmin=rmid1, rmax=rmid2, halfz=halfz_tilted)
            tiltedVols.append(self.cylVolHelper.createTrackingVolume(self.gctx, tiltedLayers2, tiltedBounds2, zshift_tilted,  f"PixelEndcap{ec}Tilted2"))

            tiltedLayers3 = [ self.layerCreator.discLayer(self.gctx, endcap[7][eta].flatten(), 1, 1) for eta in range(9) ]
            tiltedBounds3 = acts.CylinderVolumeBounds(rmin=rmid2, rmax=rmax, halfz=halfz_tilted)
            tiltedVols.append(self.cylVolHelper.createTrackingVolume(self.gctx, tiltedLayers3, tiltedBounds3, zshift_tilted,  f"PixelEndcap{ec}Tilted3"))

            pixel2Vols.append(self.cylVolHelper.createContainerTrackingVolume(self.gctx, tiltedVols))

            # Outer part of endcaps
            halfz_outer = 0.5 * (zmax_ec - zmax_tilted)
            zshift_outer = sign * (zmax_tilted + halfz_outer)

            outerVols = []

            tiltedLayers1 = [ self.layerCreator.discLayer(self.gctx, endcap[4][eta].flatten(), 1, 1) for eta in range(11) ]
            tiltedBounds1 = acts.CylinderVolumeBounds(rmin=rmin, rmax=rmid1, halfz=halfz_outer)
            outerVols.append(self.cylVolHelper.createTrackingVolume(self.gctx, tiltedLayers1, tiltedBounds1, zshift_outer,  f"PixelEndcap{ec}Outer1"))

            tiltedLayers2 = [ self.layerCreator.discLayer(self.gctx, endcap[6][eta].flatten(), 1, 1) for eta in range(8) ]
            tiltedBounds2 = acts.CylinderVolumeBounds(rmin=rmid1, rmax=rmid2, halfz=halfz_outer)
            outerVols.append(self.cylVolHelper.createTrackingVolume(self.gctx, tiltedLayers2, tiltedBounds2, zshift_outer,  f"PixelEndcap{ec}Outer2"))

            tiltedLayers3 = [ self.layerCreator.discLayer(self.gctx, endcap[8][eta].flatten(), 1, 1) for eta in range(9) ]
            tiltedBounds3 = acts.CylinderVolumeBounds(rmin=rmid2, rmax=rmax, halfz=halfz_outer)
            outerVols.append(self.cylVolHelper.createTrackingVolume(self.gctx, tiltedLayers3, tiltedBounds3, zshift_outer,  f"PixelEndcap{ec}Outer3"))

            outerPart = self.cylVolHelper.createContainerTrackingVolume(self.gctx, outerVols)
            pixel2Vols.append(self.cylVolHelper.createContainerTrackingVolume(self.gctx, outerVols))

        return self.cylVolHelper.createContainerTrackingVolume(self.gctx, pixel2Vols)


    def buildStrips(self):
        self.info("Build barrel strips")

        rmin = self.rStripInner
        rmax = self.rMax
        zmax_barrel = 1400
        zmax_ec = self.zMaxEc

        stripBarrel = self.index_hierarchy["STRIP"][0]

        stripVols = []

        stripBarrelBounds = acts.CylinderVolumeBounds(rmin=rmin, rmax=rmax, halfz=zmax_barrel)
        stripBarrelLayers = [ self.layerCreator.cylinderLayer(self.gctx, stripBarrel[lid].flatten(), 1, 1) for lid in range(4) ]
        stripVols.append(self.cylVolHelper.createTrackingVolume(self.gctx, stripBarrelLayers, stripBarrelBounds, 0.0, "StripBarrel"))

        for ec in [-2, 2]:
            self.info(f"Build strips endcaps {ec}")
            pixelEndcap = self.index_hierarchy["STRIP"][ec]

            sign = ec / abs(ec)
            halfz_ec = 0.5 * (zmax_ec - zmax_barrel)
            shift = sign * (zmax_barrel + halfz_ec)

            stripEndcapBounds = acts.CylinderVolumeBounds(rmin=rmin, rmax=rmax, halfz=halfz_ec)
            stripEndcapLayers = [ self.layerCreator.discLayer(self.gctx, pixelEndcap[lid].flatten(), 1, 1) for lid in range(6) ]
            stripVols.append(self.cylVolHelper.createTrackingVolume(self.gctx, stripEndcapLayers, stripEndcapBounds, shift, f"StripEndcap{ec}"))

        return self.cylVolHelper.createContainerTrackingVolume(self.gctx, stripVols)


    def finalize(self):
        self.info(f"Build beampipe")
        beampipeVol = acts.TrackingVolume(acts.CylinderVolumeBounds(rmin=0, rmax=self.rPixInner, halfz=self.zMaxEc), "BeamPipe")

        innerPixVol = self.buildInnerPixel()
        outerPixVol = self.buildOuterPixel()
        stripVol = self.buildStrips()

        highestVol = self.cylVolHelper.createContainerTrackingVolume(self.gctx, [beampipeVol, innerPixVol, outerPixVol, stripVol])

        trkGeoLogger = acts.logging.getLogger("TrackingGeometry")
        trkGeoLogger.setLevel(acts.logging.INFO)

        def hook(geoid, surface):
            return geoid

        self.info("Finally build tracking geometry")
        return acts.TrackingGeometry(highestVol, acts.GeometryIdentifierHook(hook), trkGeoLogger)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", type=str, default="", help="Input SQL file")
    args = p.parse_args()

    gctx = acts.GeometryContext()
    logLevel = acts.logging.INFO

    # TODO this is not yet done
    # materialDecorator = None
    # if args.map != "":
    #     print("Loading a material decorator from file:", args.map)
    #     materialDecorator = acts.IMaterialDecorator.fromFile(args.map)

    gmTree = acts.geomodel.readFromDb(args.input)

    gmFactoryConfig = gm.GeoModelDetectorSurfaceFactory.Config()
    gmFactoryConfig.shapeConverters = [
        gm.GeoBoxConverter(),
        gm.GeoTrdConverter(),
        gm.GeoIntersectionAnnulusConverter(),
    ]
    gmFactory = gm.GeoModelDetectorSurfaceFactory(gmFactoryConfig, logLevel)
    gmFactoryOptions = gm.GeoModelDetectorSurfaceFactory.Options()
    gmFactoryOptions.queries = "GeoModelXML"
    gmFactoryCache = gm.GeoModelDetectorSurfaceFactory.Cache()

    gmFactory.construct(gmFactoryCache, gctx, gmTree, gmFactoryOptions)

    itkBuilder = ItkBuilder(gmFactoryCache, gctx)

    # Uncomment this to build individual parts of the detector
    # itkBuilder.buildInnerPixel()
    # itkBuilder.buildOuterPixel()
    # itkBuilder.buildStrips()

    trkGeometry = itkBuilder.finalize()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))
    runPropagation(trkGeometry, field, outputDir="propagation").run()

if "__main__" == __name__:
    main()

