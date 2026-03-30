include(ActsTargetLinkLibrariesSystem)

acts_add_library(
    ExamplesIoCsv
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
    src/CsvMuonSpacePointReader.cpp
    src/CsvMuonSegmentReader.cpp
    src/CsvProtoTrackWriter.cpp
    src/CsvSpacePointWriter.cpp
    src/CsvGnnGraphWriter.cpp
    src/CsvGnnGraphReader.cpp
    src/CsvBFieldWriter.cpp
    src/CsvVertexWriter.cpp
)
target_include_directories(
    ActsExamplesIoCsv
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoCsv
    PRIVATE Acts::ExamplesFramework Acts::ExamplesDigitization Threads::Threads
)
