# shared code 
add_subdirectory(Common)

# Truth fitting for Generic detector 
add_executable(
  ActsExampleTruthTracksGeneric
  GenericRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksGeneric
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorGeneric) 
# Truth fitting for TGeo detector 
add_executable(
  ActsExampleTruthTracksTGeo
  TGeoRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksTGeo
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorTGeo)

# Truth fitting for Telescope detector
add_executable(
  ActsExampleTruthTracksTelescope
  TelescopeRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksTelescope
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTelescope)

# CKF reco for Generic detector 
add_executable(
  ActsExampleCKFTracksGeneric
  GenericRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksGeneric
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorGeneric)

# CKF reco for TGeo detector 
add_executable(
  ActsExampleCKFTracksTGeo
  TGeoRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksTGeo
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorTGeo)

# CKF reco for Telescope detector 
add_executable(
  ActsExampleCKFTracksTelescope
  TelescopeRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksTelescope
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTelescope)

# Hough Transform
add_executable(
  ActsExampleHough
  HoughExample.cpp)
target_link_libraries(
  ActsExampleHough
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorGeneric
)

# Chi2 fitting for Telescope detector
add_executable(
  ActsExampleChi2TracksTelescope
  TelescopeRecChi2Tracks.cpp)
target_link_libraries(
  ActsExampleChi2TracksTelescope
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTelescope)

# Seeding for Generic detector
add_executable(
  ActsExampleSeedingGeneric
  GenericSeedingExample.cpp)
target_link_libraries(
  ActsExampleSeedingGeneric
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorGeneric
)

# Seeding for TGeo detector
add_executable(
  ActsExampleSeedingTGeo
  TGeoSeedingExample.cpp)
target_link_libraries(
  ActsExampleSeedingTGeo
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTGeo
)

# =================================================
# Convert measurements to spacepoints and write sp to csv
add_executable(
  ActsExampleMeasurementsToSPTGeo
  MeasurementsToSPTGeo.cpp
)
target_link_libraries(
  ActsExampleMeasurementsToSPTGeo
  PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesIoCsv
    ActsExamplesIoJson
    ActsExamplesIoRoot
    ActsExamplesDetectorTGeo
    ActsExamplesMagneticField
    ActsExamplesTrackFinding
    ActsExamplesRecTracksCommon
)

add_executable(
  ActsExampleMeasurementsToSPGeneric
  MeasurementsToSPGeneric.cpp)
target_link_libraries(
  ActsExampleMeasurementsToSPGeneric
    PRIVATE
      ActsExamplesFramework
      ActsExamplesCommon
      ActsExamplesIoCsv
      ActsExamplesIoJson
      ActsExamplesIoRoot
      ActsExamplesDetectorGeneric
      ActsExamplesMagneticField
      ActsExamplesTrackFinding
      ActsExamplesRecTracksCommon)
# =================================================

install(
  TARGETS
    ActsExampleTruthTracksGeneric
    ActsExampleTruthTracksTGeo
    ActsExampleTruthTracksTelescope
    ActsExampleCKFTracksGeneric
    ActsExampleCKFTracksTGeo
    ActsExampleCKFTracksTelescope
    ActsExampleChi2TracksTelescope
    ActsExampleSeedingGeneric
    ActsExampleHough
    ActsExampleSeedingTGeo
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
