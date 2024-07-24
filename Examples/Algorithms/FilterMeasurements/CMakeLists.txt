add_library(
  ActsExamplesFilterMeasurements SHARED
  src/FilterMeasurementsAlgorithm.cpp
)

target_include_directories(
  ActsExamplesFilterMeasurements
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
  ActsExamplesFilterMeasurements
  PUBLIC
    ActsCore
    ActsExamplesFramework
)

install(
  TARGETS ActsExamplesFilterMeasurements
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
