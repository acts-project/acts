add_library(
  ActsExamplesSeeding INTERFACE)
target_include_directories(
  ActsExamplesSeeding
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesSeeding
  INTERFACE ActsCore ActsExamplesFramework)

# interface libraries do not exist in the filesystem; no installation needed