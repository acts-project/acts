add_executable(
  ActsExampleReadCsvGeneric
  GenericReadCsvExample.cpp)
target_link_libraries(
  ActsExampleReadCsvGeneric
  PRIVATE
    ActsCore
    ActsExamplesFramework ActsExamplesCommon ActsExamplesIoCsv
    ActsExamplesPrinters ActsExamplesDetectorGeneric
    Boost::program_options)

install(
  TARGETS ActsExampleReadCsvGeneric
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
