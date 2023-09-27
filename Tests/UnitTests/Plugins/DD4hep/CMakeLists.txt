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

set(DD4HEP_LIBRARY_PATH "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}:${DD4hep_DIR}/${CMAKE_INSTALL_LIBDIR}")

if (EXISTS ${DD4hep_DIR}/${CMAKE_INSTALL_LIBDIR})
  # set test library dependencies
  set(unittest_extra_libraries ActsTestsDD4hepFactories ActsPluginDD4hep)
  set(dd4hep_tests DD4hepDetectorElement DD4hepCylinderLayerStructure)
  foreach(_test ${dd4hep_tests})
    add_unittest(${_test} ${_test}Tests.cpp)
    if(APPLE)
      set_property(TEST DD4hepDetectorElement PROPERTY ENVIRONMENT "DYLD_LIBRARY_PATH=${DD4HEP_LIBRARY_PATH}") 
    else()
      set_property(TEST DD4hepDetectorElement PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${DD4HEP_LIBRARY_PATH}")
    endif()
  endforeach()
else()
    message(WARNING "DD4hep libraries not found under ${DD4hep_DIR}/${CMAKE_INSTALL_LIBDIR}, DD4hep tests will not be built.")
endif()