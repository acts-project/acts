add_library(
  ActsExamplesIoPerformance SHARED
  ActsExamples/Io/Performance/TrackFinderPerformanceWriter.cpp
  ActsExamples/Io/Performance/TrackFitterPerformanceWriter.cpp
  ActsExamples/Io/Performance/CKFPerformanceWriter.cpp
  ActsExamples/Io/Performance/SeedingPerformanceWriter.cpp  
  )
target_include_directories(
  ActsExamplesIoPerformance
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesIoPerformance
  PUBLIC ActsExamplesFramework
  PRIVATE ActsCore ActsExamplesSeeding  ActsPluginDigitization ROOT::Core ROOT::Tree)

install(
  TARGETS ActsExamplesIoPerformance
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
