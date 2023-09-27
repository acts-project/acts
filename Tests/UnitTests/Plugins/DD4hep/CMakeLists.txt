add_library(
    ActsTestsDD4hepFactories SHARED
    DD4hepTestsHelper.cpp
    LayerFactory.cpp
    PrimitivesFactory.cpp)

target_include_directories(
    ActsTestsDD4hepFactories
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
    ActsTestsDD4hepFactories PUBLIC ActsPluginDD4hep DD4hep::DDCore)

# Not sure why this needs to be set, but dd4hep fails to compile otherwise
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(ActsTestsDD4hepFactories PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
dd4hep_set_version(ActsTestsDD4hepFactories MAJOR 1 MINOR 0 PATCH 0)
dd4hep_generate_rootmap(ActsTestsDD4hepFactories)

install(
  TARGETS ActsTestsDD4hepFactories
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

set(DD4HEP_LIBRARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}:${DD4hep_DIR}/${CMAKE_INSTALL_LIBDIR}")

if (EXISTS ${DD4hep_DIR}/${CMAKE_INSTALL_LIBDIR})
  # set test library dependencies
  set(unittest_extra_libraries ActsTestsDD4hepFactories ActsPluginDD4hep)
  set(dd4hep_tests DD4hepDetectorElement DD4hepCylinderLayerStructure)
  foreach(_test ${dd4hep_tests})
    add_unittest(${_test} ${_test}Tests.cpp)
    if(APPLE)
      message(STATUS "Setting DYLD_LIBRARY_PATH=${DD4HEP_LIBRARY_PATH} for test ${_test}")
      set_property(TEST ${_test} PROPERTY ENVIRONMENT "DYLD_LIBRARY_PATH=${DD4HEP_LIBRARY_PATH}") 
    else()
      set_property(TEST ${_test} PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${DD4HEP_LIBRARY_PATH}")
    endif()
  endforeach()
else()
    message(WARNING "DD4hep libraries not found under ${DD4hep_DIR}/${CMAKE_INSTALL_LIBDIR}, DD4hep tests will not be built.")
endif()