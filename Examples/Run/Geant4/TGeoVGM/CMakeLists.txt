add_library(ActsExamplesTGeoVGM 
    SHARED 
    TGeoVGMDetectorConstruction.cpp)

target_compile_definitions(
  ActsExamplesTGeoVGM
  PUBLIC ${Geant4_DEFINITIONS})
target_include_directories(
  ActsExamplesTGeoVGM
  SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS} ${VGM_INCLUDE_DIRS})
target_include_directories(
  ActsExamplesTGeoVGM
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesTGeoVGM
  PUBLIC
    ActsCore ActsExamplesFramework
    ${Geant4_LIBRARIES} ${VGM_LIBRARIES} ROOT::Geom)

add_executable(
  ActsExampleGeantinoRecordingTGeoVGM
  TGeoVGMGeantinoRecording.cpp)
target_link_libraries(
  ActsExampleGeantinoRecordingTGeoVGM
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDetectorDD4hep
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesTGeoVGM
    ActsExamplesIoRoot
    Boost::program_options)

install(
  TARGETS ActsExampleGeantinoRecordingTGeoVGM
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})