message(STATUS "Building Boost ${_acts_boost_version}")

include(ExternalProject)

string(REPLACE "." "_" _acts_boost_version_ ${_acts_boost_version})

set(Boost_URL "https://boostorg.jfrog.io/artifactory/main/release/${_acts_boost_version}/source/boost_${_acts_boost_version_}.tar.gz")
set(Boost_INSTALL_DIR install)
set(Boost_INCLUDE_DIR ${Boost_INSTALL_DIR}/include)
set(Boost_LIBRARY_DIR ${Boost_INSTALL_DIR}/lib)

if(ACTS_BUILD_UNITTESTS OR ACTS_BUILD_EXAMPLES)
  set(Boost_LIBRARY_LIST "filesystem,program_options")
  set(Boost_BUILD_FILESYSTEM true)
  set(Boost_BUILD_PROGRAM_OPTIONS true)
endif()

if(ACTS_BUILD_UNITTESTS)
  string(APPEND Boost_LIBRARY_LIST ",test")
  set(Boost_BUILD_UNIT_TEST_FRAMEWORK true)
endif()

if(DEFINED Boost_LIBRARY_LIST)
  set(Boost_CONFIGURE_COMMAND ./bootstrap.sh --with-libraries=${Boost_LIBRARY_LIST} --prefix=${Boost_INSTALL_DIR})
  set(Boost_BUILD_COMMAND ./b2)
  set(Boost_INSTALL_COMMAND ./b2 install)
else()
  set(Boost_CONFIGURE_COMMAND "true")
  set(Boost_BUILD_COMMAND "true")
  set(Boost_INSTALL_COMMAND sh -c "mkdir -p ${Boost_INSTALL_DIR} && cp -r boost ${Boost_INSTALL_DIR}")
endif()
  
ExternalProject_Add(
    boost_build
    URL ${Boost_URL}
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ${Boost_CONFIGURE_COMMAND}
    BUILD_COMMAND ${Boost_BUILD_COMMAND}
    INSTALL_COMMAND ${Boost_INSTALL_COMMAND}
    LOG_CONFIGURE ON
    LOG_BUILD ON
    LOG_INSTALL ON
)

set(Boost_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/boost_build-prefix/src/boost_build/install)
set(Boost_INCLUDE_DIR ${Boost_INSTALL_DIR}/include)
set(Boost_LIBRARY_DIR ${Boost_INSTALL_DIR}/lib)

include_directories(${Boost_INCLUDE_DIR})
add_library(Boost::boost INTERFACE IMPORTED GLOBAL)
add_dependencies(Boost::boost boost_build)

if(Boost_BUILD_PROGRAM_OPTIONS)
  add_library(boost_program_options SHARED IMPORTED GLOBAL)
  set_target_properties(boost_program_options PROPERTIES IMPORTED_LOCATION ${Boost_LIBRARY_DIR}/libboost_program_options.so.${_acts_boost_version})
  add_dependencies(boost_program_options boost_build)
  add_library(Boost::program_options ALIAS boost_program_options)
  install(FILES ${Boost_LIBRARY_DIR}/libboost_program_options.so.${_acts_boost_version} DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(Boost_BUILD_FILESYSTEM)
  add_library(boost_filesystem SHARED IMPORTED GLOBAL)
  set_target_properties(boost_filesystem PROPERTIES IMPORTED_LOCATION ${Boost_LIBRARY_DIR}/libboost_filesystem.so.${_acts_boost_version})
  add_dependencies(boost_filesystem boost_build)
  add_library(Boost::filesystem ALIAS boost_filesystem)
  install(FILES ${Boost_LIBRARY_DIR}/libboost_filesystem.so.${_acts_boost_version} DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(Boost_BUILD_UNIT_TEST_FRAMEWORK)
  add_library(boost_unit_test_framework SHARED IMPORTED GLOBAL)
  set_target_properties(boost_unit_test_framework PROPERTIES IMPORTED_LOCATION ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.so.${_acts_boost_version})
  add_dependencies(boost_unit_test_framework boost_build)
  add_library(Boost::unit_test_framework ALIAS boost_unit_test_framework)
  install(FILES ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.so.${_acts_boost_version} DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

