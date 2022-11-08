add_library(
  ActsExamplesCommon SHARED
  src/CommonGeometry.cpp
  src/CommonOptions.cpp
  src/CommonSimulation.cpp
  src/GeometryExampleBase.cpp
  src/CommonMaterialMapping.cpp
  src/MaterialValidationBase.cpp
  src/PropagationExampleBase.cpp
  src/ReconstructionBase.cpp
  src/TrackFindingOptions.cpp
  src/TrackFittingChi2Options.cpp
  src/TrackFittingOptions.cpp
  src/ParticleSmearingOptions.cpp
  src/ParticleSelectorOptions.cpp
  src/ParticleGunOptions.cpp
  src/DigitizationOptions.cpp
  src/MagneticFieldOptions.cpp
  src/CsvOptionsReader.cpp
  src/CsvOptionsWriter.cpp
  src/JsonOptionsWriter.cpp
  src/VertexingOptions.cpp
  src/SpacePointMakerOptions.cpp
  src/TruthSeedSelectorOptions.cpp
  src/NuclearInteractionOptions.cpp)
target_include_directories(
  ActsExamplesCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesCommon
  PUBLIC
    ActsCore
    ActsExamplesFramework ActsExamplesMagneticField
    ActsExamplesDetectorsCommon ActsExamplesPropagation ActsExamplesGenerators
    ActsExamplesMaterialMapping ActsExamplesIoCsv ActsExamplesIoJson
    ActsExamplesIoRoot ActsExamplesIoObj ActsExamplesIoPerformance
    ActsExamplesDigitization ActsExamplesTrackFitting
    ActsExamplesTrackFittingChi2
    ActsExamplesTrackFinding ActsExamplesTruthTracking
    ActsExamplesIoCsv
  PRIVATE
    dfelibs)

install(
  TARGETS ActsExamplesCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
