import acts
import acts.examples
import acts.examples.detector as detector
from acts import LayerStructureBuilder
import math

# main script runs
jsOptions = acts.examples.SurfaceJsonOptions()
jsOptions.inputFile = "odd-input.json"

# Where to pick the surfaces from
surfacesHierarchyMap = acts.examples.readSurfaceFromJson(jsOptions)
svlMap = acts.examples.extractVolumeLayerSurfaces(surfacesHierarchyMap, True)

LayerSurfaces = acts.LayerStructureBuilder.SurfacesHolder

# Create the detector structure
detectorRmin = 0.0
detectorRmax = 1200
detectorZmin = -3100
detectorZmax = -detectorZmin
beamPipeRmax = 27.0

nominal_transform = acts.Transform3([0.0, 0.0, 0.0])

# Beam pipe section ###########################################################
beamPipe = detector.DetectorVolume(
    "BeamPipe",
    acts.VolumeBoundsType.Cylinder,
    [detectorRmin, beamPipeRmax, detectorZmax],
    nominal_transform,
)

# Pixel section ###############################################################
pixelRmax = 200.0
pixelZmid = 580

# Endcap configurations
pixelEndcapConfigs = [
    [620.0, 5.0],
    [720.0, 5.0],
    [840.0, 5.0],
    [980.0, 5.0],
    [1120.0, 5.0],
    [1320.0, 5.0],
    [1520.0, 5.0],
]
nPEs = len(pixelEndcapConfigs)

# Pixel endcap section
pixelEndcapBinningR = acts.ProtoBinning(
    acts.Binning.r, acts.Binning.bound, 40, 175, 2, 1
)
pixelEndcapBinningPhi = acts.ProtoBinning(
    acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 56, 1
)
pixelEndcapConfigsN = [
    [-pixelEndcapConfigs[i][0], pixelEndcapConfigs[i][1]] for i in range(nPEs)
]
pixelEndcapConfigsN.reverse()

pixelNegativeEndcap = detector.createEndcap(
    "PixelNegativeEndcap",
    [beamPipeRmax, pixelRmax, -detectorZmax, -pixelZmid],
    pixelEndcapConfigsN,
    [LayerSurfaces(svlMap[16][4 + i * 2]) for i in range(nPEs)],
    [[pixelEndcapBinningR, pixelEndcapBinningPhi] for i in range(nPEs)],
)

pixelPositiveEndcap = detector.createEndcap(
    "PixelPositiveEndcap",
    [beamPipeRmax, pixelRmax, pixelZmid, detectorZmax],
    pixelEndcapConfigs,
    [LayerSurfaces(svlMap[18][2 + i * 2]) for i in range(nPEs)],
    [[pixelEndcapBinningR, pixelEndcapBinningPhi] for i in range(nPEs)],
)

# Pixel barrel section
pixelBarrelDims = [beamPipeRmax, pixelRmax, -pixelZmid, pixelZmid]
nPBs = len(pixelBarrelDims)

pixelBarrelConfigs = [[34.0, 5.0], [70.0, 5], [116.0, 5.0], [172.0, 5.0]]
pixelBarrelSurfaces = [LayerSurfaces(svlMap[17][2 + i * 2]) for i in range(nPBs)]
pixelBarrelBinningZ = acts.ProtoBinning(
    acts.Binning.z, acts.Binning.bound, -500, 500, 14, 1
)
pixelBarrelBinnings = [
    [pixelBarrelBinningZ, detector.phiBinning(16, 1)],
    [pixelBarrelBinningZ, detector.phiBinning(32, 1)],
    [pixelBarrelBinningZ, detector.phiBinning(52, 1)],
    [pixelBarrelBinningZ, detector.phiBinning(78, 1)],
]

pixelBarrel = detector.createBarrel(
    "PixelBarrel",
    pixelBarrelDims,
    pixelBarrelConfigs,
    pixelBarrelSurfaces,
    pixelBarrelBinnings,
)

pixelContainer = detector.ContainerStructure(
    "Pixel", [pixelNegativeEndcap, pixelBarrel, pixelPositiveEndcap], acts.Binning.z
)

# PST #########################################################################
pstRmax = 220.0
pst = detector.DetectorVolume(
    "PST",
    acts.VolumeBoundsType.Cylinder,
    [pixelRmax, pstRmax, detectorZmax],
    nominal_transform,
)

# Short strip section #########################################################
sstripRmax = 720.0
sstripZmid = 1250.0

# SStrip negative/positive endcap
sstripEndcapBinningR = acts.ProtoBinning(
    acts.Binning.r, acts.Binning.bound, 230, 710, 3, 1
)
sstripEndcapBinningPhi = acts.ProtoBinning(
    acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 40, 1
)
sstripEndcapConfigs = [
    [1300.0, 15.0],
    [1550.0, 15.0],
    [1850.0, 15.0],
    [2200, 15.0],
    [2550.0, 15.0],
    [2950.0, 15.0],
]
nSEs = len(sstripEndcapConfigs)

sstripEndcapConfigsN = [
    [-sstripEndcapConfigs[i][0], sstripEndcapConfigs[i][1]] for i in range(nSEs)
]
sstripEndcapConfigsN.reverse()

sstripNegativeEndcap = detector.createEndcap(
    "ShortsStripsNegativeEndcap",
    [pstRmax, sstripRmax, -detectorZmax, -sstripZmid],
    sstripEndcapConfigsN,
    [LayerSurfaces(svlMap[23][2 + i * 2]) for i in range(nSEs)],
    [[sstripEndcapBinningR, sstripEndcapBinningPhi] for i in range(nSEs)],
)

sstripPositiveEndcap = detector.createEndcap(
    "ShortStripsPositiveEndcap",
    [pstRmax, sstripRmax, sstripZmid, detectorZmax],
    sstripEndcapConfigs,
    [LayerSurfaces(svlMap[25][2 + i * 2]) for i in range(nSEs)],
    [[sstripEndcapBinningR, sstripEndcapBinningPhi] for i in range(nSEs)],
)

# SStrip barrel
sstripBarrelConfigs = [[260.0, 25], [360.0, 25.0], [500.0, 25.0], [660.0, 25.0]]
nSBs = len(sstripBarrelConfigs)

sstripBarrelBinningZ = acts.ProtoBinning(
    acts.Binning.z, acts.Binning.bound, -1100, 1100, 21, 1
)
sstripBarrelPhiBins = [40, 56, 78, 102]
sstripBarrel = detector.createBarrel(
    "ShortStripsBarrel",
    [pstRmax, sstripRmax, -sstripZmid, sstripZmid],
    sstripBarrelConfigs,
    [LayerSurfaces(svlMap[24][2 + i * 2]) for i in range(nSBs)],
    [
        [sstripBarrelBinningZ, detector.phiBinning(sstripBarrelPhiBins[i], 1)]
        for i in range(nSBs)
    ],
)

sstripContainer = detector.ContainerStructure(
    "ShortStrips",
    [sstripNegativeEndcap, sstripBarrel, sstripPositiveEndcap],
    acts.Binning.z,
)

# Long strip section ##########################################################

# LStrip section
lstripRmax = 1100.0
lstripZmid = sstripZmid

# LStrip negative/positive endcap
lstripEndcapBinningR = acts.ProtoBinning(
    acts.Binning.r, acts.Binning.bound, 720, 1020, 2, 1
)
lstripEndcapBinningPhi = acts.ProtoBinning(
    acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 40, 1
)
lstripEndcapConfigs = [
    [1300.0, 25.0],
    [1600.0, 25.0],
    [1900.0, 25.0],
    [2250.0, 25.0],
    [2600.0, 25.0],
    [3000.0, 25.0],
]
nLEs = len(lstripEndcapConfigs)

lstripEndcapConfigsN = [
    [-lstripEndcapConfigs[i][0], lstripEndcapConfigs[i][1]] for i in range(nLEs)
]
lstripEndcapConfigsN.reverse()

lstripNegativeEndcap = detector.createEndcap(
    "LongStripsNegativeEndcap",
    [sstripRmax, lstripRmax, -detectorZmax, -lstripZmid],
    lstripEndcapConfigsN,
    [LayerSurfaces(svlMap[28][2 + i * 2]) for i in range(nLEs)],
    [[lstripEndcapBinningR, lstripEndcapBinningPhi] for i in range(nLEs)],
)

lstripPositiveEndcap = detector.createEndcap(
    "LongStripsPositiveEndcap",
    [sstripRmax, lstripRmax, lstripZmid, detectorZmax],
    lstripEndcapConfigs,
    [LayerSurfaces(svlMap[30][2 + i * 2]) for i in range(nLEs)],
    [[lstripEndcapBinningR, lstripEndcapBinningPhi] for i in range(nLEs)],
)

# SStrip barrel
lstripBarrelConfigs = [[830.0, 35.0], [1030, 35.0]]
lstripBarrelBinningZ = acts.ProtoBinning(
    acts.Binning.z, acts.Binning.bound, -1100, 1100, 21, 1
)
lstripBarrelPhiBins = [40, 56]
nLBs = len(lstripBarrelConfigs)

lstripBarrel = detector.createBarrel(
    "LongStripsBarrel",
    [sstripRmax, lstripRmax, -lstripZmid, lstripZmid],
    lstripBarrelConfigs,
    [LayerSurfaces(svlMap[29][2 + i * 2]) for i in range(nLBs)],
    [
        [lstripBarrelBinningZ, detector.phiBinning(lstripBarrelPhiBins[i], 1)]
        for i in range(2)
    ],
)

lstripContainer = detector.ContainerStructure(
    "LongStrips",
    [lstripNegativeEndcap, lstripBarrel, lstripPositiveEndcap],
    acts.Binning.z,
)

### Detector section ##########################################################
detectorContainer = detector.ContainerStructure(
    "Detector",
    [beamPipe, pixelContainer, pst, sstripContainer, lstripContainer],
    acts.Binning.r,
)

detectorBuilderConf = acts.DetectorBuilder.Config()
detectorBuilderConf.name = detectorContainer.getName()
detectorBuilderConf.builder = detectorContainer.getBuilder()

detectorBuilder = acts.DetectorBuilder(
    detectorBuilderConf, "ODD Detector Builder", acts.logging.VERBOSE
)

detector = detectorBuilder.construct(acts.GeometryContext())
