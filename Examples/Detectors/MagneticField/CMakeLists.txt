add_library(
  ActsExamplesMagneticField SHARED
  src/BFieldOptions.cpp
  src/BFieldScalor.cpp
  src/BFieldUtils.cpp)
target_include_directories(
  ActsExamplesMagneticField
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesMagneticField
  PUBLIC ActsCore ActsExamplesFramework Boost::program_options ROOT::Core ROOT::Tree)

install(
  TARGETS ActsExamplesMagneticField
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
