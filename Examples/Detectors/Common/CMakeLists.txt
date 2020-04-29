add_library(
  ActsExamplesDetectorsCommon INTERFACE)
target_include_directories(
  ActsExamplesDetectorsCommon
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDetectorsCommon INTERFACE Boost::program_options)

# interface libraries do not exist in the filesystem; no installation needed
