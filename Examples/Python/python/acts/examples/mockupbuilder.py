import acts
import acts.examples.geant4 as actsG4

mockupConfig = actsG4.MockupSectorBuilder.Config()

mockupChamberConfigInner = actsG4.MockupSectorBuilder.ChamberConfig()
mockupChamberConfigInner.name = "Inner_Detector_Chamber"
mockupChamberConfigInner.SensitiveNames = ["Inner_Skin"]
mockupChamberConfigInner.PassiveNames = ["xx"]

mockupChamberConfigMiddle = actsG4.MockupSectorBuilder.ChamberConfig()
mockupChamberConfigMiddle.name = "Middle_Detector_Chamber"
mockupChamberConfigMiddle.SensitiveNames = ["Middle_Skin"]
mockupChamberConfigMiddle.PassiveNames = ["xx"]

mockupChamberConfigOuter = actsG4.MockupSectorBuilder.ChamberConfig()
mockupChamberConfigOuter.name = "Outer_Detector_Chamber"
mockupChamberConfigOuter.SensitiveNames = ["Outer_Skin"]
mockupChamberConfigOuter.PassiveNames = ["xx"]

mockupConfig.gdmlPath = "/home/dimitra/ACTS_Dimitra/acts/Examples/Python/python/acts/examples/MuonChamber.gdml"
mockupConfig.NumberOfSectors = 8

mockupBuilder = actsG4.MockupSectorBuilder(mockupConfig)
print("hello")
detectorVolumeInner = mockupBuilder.buildChamber(mockupChamberConfigInner)
print("inner built")
detectorVolumeOuter = mockupBuilder.buildChamber(mockupChamberConfigOuter)
print("outer built")
detectorVolumeMiddle = mockupBuilder.buildChamber(mockupChamberConfigMiddle)
print("middle built")
print("hello")

detectorVolumes = [detectorVolumeInner, detectorVolumeMiddle, detectorVolumeOuter]

detectorVolumeSector = mockupBuilder.buildSector(detectorVolumes)

mockupBuilder.drawSector(detectorVolumeSector, "sector_drawn_from_python")
