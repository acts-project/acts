add_library(
  ActsExamplesDigitization SHARED
  src/DigitizationAlgorithm.cpp
  src/HitSmearing.cpp)
target_include_directories(
  ActsExamplesDigitization
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDigitization
  PRIVATE
    ActsCore ActsDigitizationPlugin ActsIdentificationPlugin
    ActsExamplesFramework
    Boost::program_options)

install(
  TARGETS ActsExamplesDigitization
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
