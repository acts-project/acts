set(SOURCES 
  src/AmbiguityResolutionML.cpp
  src/AmbiguityResolutionMLAlgorithm.cpp
)

if(ACTS_BUILD_PLUGIN_MLPACK)
  list(APPEND SOURCES 
    src/AmbiguityResolutionMLDBScanAlgorithm.cpp   
  )
endif()

add_library(
  ActsExamplesTrackFindingML SHARED
  ${SOURCES}
)

target_include_directories(
  ActsExamplesTrackFindingML
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
  ActsExamplesTrackFindingML
  PUBLIC
    ActsCore
    ActsPluginOnnx
    ActsExamplesFramework
)

if(ACTS_BUILD_PLUGIN_MLPACK)
  target_link_libraries(
    ActsExamplesTrackFindingML PUBLIC ActsPluginmlpack)
endif()

install(
  TARGETS ActsExamplesTrackFindingML
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
