add_library(
  ActsExamplesUtility SHARED
  src/AliasAlgorithm.cpp)
target_include_directories(
  ActsExamplesUtility
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesUtility
  PUBLIC
  ActsCore ActsExamplesFramework Boost::program_options)

install(
  TARGETS ActsExamplesUtility
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
