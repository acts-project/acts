add_library(
	ActsExamplesSeedingCommon SHARED
	SeedingExample.cpp)

target_include_directories(
  ActsExamplesSeedingCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(
  ActsExamplesSeedingCommon
  PUBLIC
    ActsCore
    ActsExamplesFramework
    ActsExamplesPropagation
    ActsExamplesSeeding
    ActsExamplesCommon
    ActsExamplesDigitization  
    ActsExamplesMagneticField
    ActsExamplesIoObj
    ActsExamplesIoJson
    ActsExamplesIoRoot
    Boost::program_options)
