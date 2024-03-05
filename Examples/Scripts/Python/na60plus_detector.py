# import the necessary modules

import acts, acts.examples
from acts.examples import geant4 as acts_g4


from acts import (
    svg,
    logging,
    Binning,
    Extent,
    GeometryContext,
    ProtoBinning,
    LayerStructureBuilder,
    VolumeBoundsType,
    VolumeStructureBuilder,
    DetectorVolumeBuilder,
    GeometryIdGenerator,
    CuboidalContainerBuilder,
    DetectorBuilder,
    Transform3,
    Range1D,
    KdtSurfaces1D,
    KdtSurfacesProvider1D)

geoContext = acts.GeometryContext()

sensitive_matches = [ 'PixSensor' ]
passive_matches = [ 'PixStnFrame' ]

[ elements, ssurfaces, psurfaces ] = acts_g4.convertSurfaces('na60VT.gdml',  sensitive_matches, passive_matches)

logLevel = logging.VERBOSE

# Write them to an obj file
drawContext = acts.GeometryContext()
sensitiveRgb = [ 0, 150, 150 ]
passiveRgb = [ 150, 150, 0]

kdtSurfaces = KdtSurfaces1D(geoContext, ssurfaces, [Binning.z] )

layers = [   [-386, -376], [-376, -256], [-256, -246], [-246, -205], [-205, -195], [-195, -155 ], [-155, -145], [-145, -75], [-75., -65], [-65, 25 ] ]

volBuilders = []
gaps = 0
layerId = 0

for il, lrange in enumerate(layers) :
    
    surfaces = kdtSurfaces.surfaces(Range1D(lrange))

    layerStr = "Layer_"+str(4-layerId)
    if len(surfaces) == 0 :
            layerStr = "Gap_"+str(gaps)
            gaps += 1
    else :
            layerId += 1

    print ('-> Building Volume: ', layerStr, ' with surfaces: ', len(surfaces))

     # Set up the shape builder: external builder
    shapeConfig = VolumeStructureBuilder.Config()
    shapeConfig.boundsType = VolumeBoundsType.Cuboid
    shapeConfig.boundValues = [ 350., 350., 0.5 * (lrange[1] - lrange[0]) ]
    shapeConfig.transform = Transform3([0, 0, 0.5 * (lrange[1] + lrange[0])])
    shapeConfig.auxiliary = "Shape[" + layerStr + "]"
    
    layerStructure = None
    if (len(surfaces) > 0) :
        surfaceProvider = KdtSurfacesProvider1D(kdtSurfaces, Extent([[Binning.z, lrange]]))

        layerConfig = LayerStructureBuilder.Config()
        layerConfig.surfacesProvider = surfaceProvider
        layerConfig.binnings = []
        layerConfig.supports = []
        layerConfig.auxiliary = layerStr

        layerStructure = LayerStructureBuilder(
                layerConfig, layerConfig.auxiliary, logLevel)
    
    volConfig = acts.DetectorVolumeBuilder.Config()
    volConfig.name = "Volume["+layerStr+"]"
    volConfig.auxiliary = volConfig.name
    volConfig.externalsBuilder = VolumeStructureBuilder(
            shapeConfig, shapeConfig.auxiliary, logLevel)    
    volConfig.internalsBuilder = layerStructure
    
    volBuilder = DetectorVolumeBuilder(volConfig, volConfig.name, logLevel)
    volBuilders += [volBuilder]

print ('-> Building Detector from ', len(volBuilders), ' volumes')

geoIdConfig = GeometryIdGenerator.Config()
geoIdConfig.containerMode = True
geoId = GeometryIdGenerator(geoIdConfig, 'GeometryIdGenerator', logLevel)

ccConfig = CuboidalContainerBuilder.Config()
ccConfig.builders = volBuilders
ccConfig.binning = Binning.z
ccConfig.geoIdGenerator = geoId
ccConfig.auxiliary = "SiliconTracker"

ccBuilder = CuboidalContainerBuilder(ccConfig, ccConfig.auxiliary, logLevel)


detConfig = DetectorBuilder.Config()
detConfig.name = "na60+"
detConfig.builder = ccBuilder
detConfig.auxiliary = detConfig.name
detConfig.geoIdGenerator = geoId

detBuilder = DetectorBuilder(detConfig, detConfig.name, logLevel)
detector = detBuilder.construct(geoContext)




