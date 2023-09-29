add_library(
    ActsExamplesDD4hepDetectorWithOptions SHARED
    src/DD4hepDetectorWithOptions.cpp
)

target_link_libraries(
    ActsExamplesDD4hepDetectorWithOptions PUBLIC
    ActsExamplesCommon ActsExamplesDetectorDD4hep)

target_include_directories(
    ActsExamplesDD4hepDetectorWithOptions
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

install(
  TARGETS ActsExamplesDD4hepDetectorWithOptions
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
