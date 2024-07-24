add_library(
  ActsExamplesMergeContainers SHARED
  src/MergeContainersAlgorithm.cpp
)

target_include_directories(
  ActsExamplesMergeContainers
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
  ActsExamplesMergeContainers
  PUBLIC
    ActsCore
    ActsExamplesFramework
)

install(
  TARGETS ActsExamplesMergeContainers
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
