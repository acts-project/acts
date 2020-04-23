add_library(
  ActsExamplesPropagation INTERFACE)
target_include_directories(
  ActsExamplesPropagation
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesPropagation
  INTERFACE ActsCore ActsExamplesFramework)

# interface libraries do not exist in the filesystem; no installation needed
