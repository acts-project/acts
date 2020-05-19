add_library(
  ActsExamplesIoPerformance SHARED
  ACTFW/Io/Performance/TrackFinderPerformanceWriter.cpp
  ACTFW/Io/Performance/TrackFitterPerformanceWriter.cpp
  ACTFW/Io/Performance/CKFPerformanceWriter.cpp)
target_include_directories(
  ActsExamplesIoPerformance
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesIoPerformance
  PUBLIC ActsExamplesFramework
  PRIVATE ActsCore ROOT::Core ROOT::Tree)

install(
  TARGETS ActsExamplesIoPerformance
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
