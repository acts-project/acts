add_library(
  ActsExamplesCalibration
  SHARED
  src/MeasurementCalibration.cpp)

target_link_libraries(
  ActsExamplesCalibration
  PUBLIC ActsCore ActsExamplesFramework Eigen3::Eigen)

target_include_directories(
  ActsExamplesCalibration
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

install(
  TARGETS ActsExamplesCalibration
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

