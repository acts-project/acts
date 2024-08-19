include(ActsTargetLinkLibrariesSystem)

add_library(
  ActsExamplesIoCsv SHARED
  src/CsvMeasurementReader.cpp
  src/CsvMeasurementWriter.cpp
  src/CsvParticleReader.cpp
  src/CsvParticleWriter.cpp
  src/CsvSimHitReader.cpp
  src/CsvSimHitWriter.cpp
  src/CsvSpacePointReader.cpp
  src/CsvTrackingGeometryWriter.cpp
  src/CsvTrackParameterReader.cpp
  src/CsvTrackParameterWriter.cpp
  src/CsvSeedWriter.cpp
  src/CsvTrackWriter.cpp
  src/CsvDriftCircleReader.cpp
  src/CsvMuonSimHitReader.cpp
  src/CsvProtoTrackWriter.cpp
  src/CsvSpacePointWriter.cpp
  src/CsvExaTrkXGraphWriter.cpp
  src/CsvExaTrkXGraphReader.cpp
  src/CsvBFieldWriter.cpp)
target_include_directories(
  ActsExamplesIoCsv
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoCsv
  PRIVATE
    ActsCore
    ActsExamplesFramework ActsExamplesDigitization
    Threads::Threads)

if(ACTS_BUILD_PLUGIN_GEOMODEL)
target_link_libraries(ActsExamplesIoCsv PRIVATE ActsPluginGeoModel)
target_compile_definitions(ActsExamplesIoCsv PRIVATE HAVE_GEOMODEL)
endif()

install(
  TARGETS ActsExamplesIoCsv
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
