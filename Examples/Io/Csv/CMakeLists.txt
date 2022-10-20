add_library(
  ActsExamplesIoCsv SHARED
  src/CsvMeasurementReader.cpp
  src/CsvMeasurementWriter.cpp
  src/CsvParticleReader.cpp
  src/CsvParticleWriter.cpp
  src/CsvPlanarClusterReader.cpp
  src/CsvPlanarClusterWriter.cpp
  src/CsvSimHitReader.cpp
  src/CsvSimHitWriter.cpp
  src/CsvSpacePointReader.cpp
  src/CsvTrackingGeometryWriter.cpp
  src/CsvMultiTrajectoryWriter.cpp 
  src/CsvSpacePointWriter.cpp)
target_include_directories(
  ActsExamplesIoCsv
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoCsv
  PRIVATE
    ActsCore ActsPluginIdentification
    ActsExamplesFramework ActsExamplesDigitization
    Threads::Threads dfelibs)

if( ACTS_BUILD_EXAMPLES_BINARIES )
  target_link_libraries(
    ActsExamplesIoCsv
    PRIVATE
    Boost::program_options
  )
endif()

install(
  TARGETS ActsExamplesIoCsv
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
