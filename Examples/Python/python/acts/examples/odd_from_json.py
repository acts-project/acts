import acts
import acts.examples
import acts.examples.detector as detector
from acts import LayerStructureBuilder
import math

# main script runs
jsOptions = acts.examples.SurfaceJsonOptions()
jsOptions.inputFile = 'odd-input.json'

# Where to pick the surfaces from
surfacesHierarchyMap = acts.examples.readSurfaceFromJson(jsOptions)
svlMap = acts.examples.extractVolumeLayerSurfaces(surfacesHierarchyMap, True)

LayerSurfaces = acts.LayerStructureBuilder.SurfacesHolder

# Create the detector structure
detectorRmin = 0.
detectorRmax = 1200
detectorZmin = -3100
detectorZmax = -detectorZmin
beamPipeRmax = 27.

nominal_transform = acts.Transform3([0., 0., 0.])

# Beam pipe section ###########################################################
beamPipe = detector.DetectorVolume(
    'BeamPipe', 
    acts.VolumeBoundsType.Cylinder, 
    [detectorRmin, beamPipeRmax, detectorZmax], 
    nominal_transform)

# Pixel section ###############################################################
pixelRmax = 200.
pixelZmid = 580

# Endcap configurations
pixelEndcapConfigs = [ [620., 5.], [720., 5.], [840., 5.], [980., 5.], [1120., 5.], [1320., 5.], [1520., 5.]]

# Pixel endcap section
pixelEndcapBinningR = acts.ProtoBinning(
    acts.Binning.r, acts.Binning.bound, 40, 175, 2, 1)
pixelEndcapBinningPhi = acts.ProtoBinning(
    acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 56, 1)
pixelEndcapConfigsN = [ [-pixelEndcapConfigs[i][0], pixelEndcapConfigs[i][1]] for i in range(7) ]
pixelEndcapConfigsN.reverse()

pixelNegativeEndcap = detector.createEndcap( 'PixelNegativeEndcap', 
                                     [ beamPipeRmax, pixelRmax, -detectorZmax, -pixelZmid ],  
                                     pixelEndcapConfigsN,
                                     [ LayerSurfaces(svlMap[16][4 + i *2]) for i in range(7)], 
                                     [ [pixelEndcapBinningR, pixelEndcapBinningPhi] for i in range (7)])

pixelPositiveEndcap = detector.createEndcap( 'PixelPositiveEndcap', 
                                     [ beamPipeRmax, pixelRmax, pixelZmid, detectorZmax ],  
                                     pixelEndcapConfigs,
                                      [ LayerSurfaces(svlMap[18][2+ i *2]) for i in range(7)], 
                                     [ [pixelEndcapBinningR, pixelEndcapBinningPhi] for i in range (7)])

# Pixel barrel section
pixelBarrelDims = [ beamPipeRmax, pixelRmax, -pixelZmid, pixelZmid]
pixelBarrelConfigs  = [ [34., 5.], [70., 5], [116., 5.], [172., 5.]]
pixelBarrelSurfaces = [ LayerSurfaces(svlMap[17][2 + i *2]) for i in range(4)]
pixelBarrelBinningZ = acts.ProtoBinning(
    acts.Binning.z, acts.Binning.bound, -500, 500, 14, 1)
pixelBarrelBinnings = [ [pixelBarrelBinningZ, detector.phiBinning(16,1)],
                      [pixelBarrelBinningZ, detector.phiBinning(32,1)], 
                      [pixelBarrelBinningZ, detector.phiBinning(52,1)], 
                      [pixelBarrelBinningZ, detector.phiBinning(78,1)]]

pixelBarrel = detector.createBarrel('PixelBarrel', 
                                    pixelBarrelDims, 
                                    pixelBarrelConfigs, 
                                    pixelBarrelSurfaces,
                                    pixelBarrelBinnings)

pixelContainer = detector.ContainerStructure('Pixel', 
                                    [ pixelNegativeEndcap, pixelBarrel, pixelPositiveEndcap ], 
                                    acts.Binning.z)

# PST #########################################################################
pstRmax = 220.
pst = detector.DetectorVolume(
        'PST', 
        acts.VolumeBoundsType.Cylinder, 
        [pixelRmax, pstRmax, detectorZmax], 
        nominal_transform)

# Short strip section #########################################################
sstripRmax = 720.
sstripZmid = 1250.

# SStrip negative/positive endcap
sstripEndcapBinningR = acts.ProtoBinning(
    acts.Binning.r, acts.Binning.bound, 230, 710, 3, 1)
sstripEndcapBinningPhi = acts.ProtoBinning(
    acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 40, 1)
sstripEndcapConfigs = [ [1300., 15.], 
                       [1550., 15.], 
                       [1850., 15.], 
                       [2200, 15.], 
                       [2550., 15.], 
                       [2950., 15.] ]

sstripEndcapConfigsN = [ [-sstripEndcapConfigs[i][0], sstripEndcapConfigs[i][1]] for i in range(6) ]
sstripEndcapConfigsN.reverse()

sstripNegativeEndcap = detector.createEndcap( 'ShortsStripsNegativeEndcap', 
                                     [ pstRmax, sstripRmax, -detectorZmax, -sstripZmid ],  
                                     sstripEndcapConfigsN,
                                     [ LayerSurfaces(svlMap[23][2 + i *2]) for i in range(6)], 
                                     [ [sstripEndcapBinningR, sstripEndcapBinningPhi] for i in range (6)])

sstripPositiveEndcap = detector.createEndcap( 'ShortStripsPositiveEndcap', 
                                     [ pstRmax, sstripRmax, sstripZmid, detectorZmax ],  
                                     sstripEndcapConfigs,
                                      [ LayerSurfaces(svlMap[25][2+ i *2]) for i in range(6)], 
                                     [ [sstripEndcapBinningR, sstripEndcapBinningPhi] for i in range (6)])

# SStrip barrel
sstripBarrelConfigs = [ [260., 25], [360., 25.], [500., 25.], [660., 25.] ]
sstripBarrelBinningZ = acts.ProtoBinning(
    acts.Binning.z, acts.Binning.bound, -1100, 1100, 21, 1)
sstripBarrelPhiBins = [40, 56, 78, 102]
sstripBarrel = detector.createBarrel('ShortStripsBarrel', 
                            [pstRmax, sstripRmax, -sstripZmid, sstripZmid],
                            sstripBarrelConfigs,
                            [ LayerSurfaces(svlMap[24][2 + i *2]) for i in range(4) ],
                            [ [sstripBarrelBinningZ, detector.phiBinning(sstripBarrelPhiBins[i],1)] for i in range (4)])


sstripContainer = detector.ContainerStructure('ShortStrips',
                                    [sstripNegativeEndcap, sstripBarrel, sstripPositiveEndcap],
                                    acts.Binning.z)

# Long strip section ##########################################################

# LStrip section
lstripRmax = 1100.
lstripZmid = sstripZmid

# LStrip negative/positive endcap
lstripEndcapBinningR = acts.ProtoBinning(
    acts.Binning.r, acts.Binning.bound, 720, 1020, 2, 1)
lstripEndcapBinningPhi = acts.ProtoBinning(
    acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 40, 1)
lstripEndcapConfigs= [ [1300., 25.], 
                       [1600., 25.], 
                       [1900., 25.], 
                       [2250., 25.], 
                       [2600., 25.], 
                       [3000., 25.]]

lstripEndcapConfigsN = [ [-lstripEndcapConfigs[i][0], lstripEndcapConfigs[i][1]] for i in range(6) ]
lstripEndcapConfigsN.reverse()

lstripNegativeEndcap = detector.createEndcap( 'LongStripsNegativeEndcap', 
                                     [ sstripRmax, lstripRmax, -detectorZmax, -lstripZmid ],  
                                     lstripEndcapConfigsN,
                                     [ LayerSurfaces(svlMap[28][2 + i *2]) for i in range(6)], 
                                     [ [lstripEndcapBinningR, lstripEndcapBinningPhi] for i in range (6)])

lstripPositiveEndcap = detector.createEndcap( 'LongStripsPositiveEndcap', 
                                     [ sstripRmax, lstripRmax, lstripZmid, detectorZmax ],  
                                     lstripEndcapConfigs,
                                      [ LayerSurfaces(svlMap[30][2+ i *2]) for i in range(6)], 
                                     [ [lstripEndcapBinningR, lstripEndcapBinningPhi] for i in range (6)])

# SStrip barrel
lstripBarrelConfigs = [ [830., 35.], [1030, 35.] ]
lstripBarrelBinningZ = acts.ProtoBinning(
    acts.Binning.z, acts.Binning.bound, -1100, 1100, 21, 1)
lstripBarrelPhiBins = [40, 56]

lstripBarrel = detector.createBarrel('LongStripsBarrel', 
                            [sstripRmax, lstripRmax, -lstripZmid, lstripZmid],
                            lstripBarrelConfigs,
                            [ LayerSurfaces(svlMap[29][2 + i *2]) for i in range(2) ],
                            [ [lstripBarrelBinningZ, detector.phiBinning(lstripBarrelPhiBins[i],1)] for i in range (2)])

lstripContainer = detector.ContainerStructure('LongStrips',
                                    [lstripNegativeEndcap, lstripBarrel, lstripPositiveEndcap],
                                    acts.Binning.z)

### Detector section ##########################################################
detectorContainer = detector.ContainerStructure('Detector', 
                                                [ beamPipe, pixelContainer, pst, sstripContainer, lstripContainer ],
                                                acts.Binning.r)

detectorBuilderConf = acts.DetectorBuilder.Config()
detectorBuilderConf.name = detectorContainer.getName()
detectorBuilderConf.builder = detectorContainer.getBuilder()

detectorBuilder = acts.DetectorBuilder(
    detectorBuilderConf, 'ODD Detector Builder', acts.logging.VERBOSE)

detector = detectorBuilder.construct(acts.GeometryContext())
