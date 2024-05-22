message(STATUS "Building Boost ${_acts_boost_recommended_version}")

include(ExternalProject)

string(REPLACE "." "_" _acts_boost_recommended_version_ ${_acts_boost_recommended_version})
set(Boost_LIBRARY_LIST "filesystem,program_options,test")
set(Boost_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/boost-prefix/src/boost)
set(Boost_INCLUDE_DIR ${Boost_SOURCE_DIR})
set(Boost_LIBRARY_DIR ${Boost_SOURCE_DIR}/stage/lib)
set(Boost_B2_CXXFLAGS -std=c++17\ -fPIC\ -O3)
set(Boost_B2_FLAGS variant=release threading=multi runtime-link=shared link=shared cxxflags=${Boost_B2_CXXFLAGS})

# Target for boost headers

# Note: Since all builds below depend on the headers, use this step to
# download the boost release

ExternalProject_Add( boost "${ACTS_BOOST_SOURCE};CONFIGURE_COMMAND;;BUILD_COMMAND;;INSTALL_COMMAND;" )

add_library(Boost::boost INTERFACE IMPORTED GLOBAL)
add_dependencies(Boost::boost boost)
target_include_directories(Boost::boost INTERFACE ${Boost_INCLUDE_DIR})
install(DIRECTORY ${Boost_INCLUDE_DIR}/boost DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Target for bootstrap step

# Note: Required by all builds but not by header-only target, so need
# a separate target

# Note: Re-use ExternalProject_Add instead of Add_Step because I can't
# figure out how to avoid race conditions and duplicate downloads with
# Add_Step targets.  From cmake's point of view they will all be
# different project but they all use the same working directory, so it
# does the right thing.  (L-G)

ExternalProject_Add(
  boost_bootstrap
  DOWNLOAD_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ./bootstrap.sh --with-libraries=${Boost_LIBRARY_LIST}
  INSTALL_COMMAND ""
  DEPENDS boost
  LOG_BUILD true
  SOURCE_DIR ${Boost_SOURCE_DIR}
  BUILD_IN_SOURCE true
  EXCLUDE_FROM_ALL true
)

# Target for Boost::program_options

ExternalProject_Add(
  boost_build_program_options
  DOWNLOAD_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ./b2 ${Boost_B2_FLAGS} --with-program_options
  INSTALL_COMMAND ""
  DEPENDS boost_bootstrap
  SOURCE_DIR ${Boost_SOURCE_DIR}
  LOG_BUILD true
  BUILD_IN_SOURCE true
  BUILD_BYPRODUCTS ${Boost_LIBRARY_DIR}/libboost_program_options.so.${_acts_boost_recommended_version}
  EXCLUDE_FROM_ALL true
)

add_library(boost_program_options SHARED IMPORTED GLOBAL)
target_link_libraries(boost_program_options INTERFACE Boost::boost)
set_target_properties(boost_program_options PROPERTIES IMPORTED_LOCATION ${Boost_LIBRARY_DIR}/libboost_program_options.so.${_acts_boost_recommended_version})
add_dependencies(boost_program_options boost_build_program_options)
add_library(Boost::program_options ALIAS boost_program_options)
install(FILES ${Boost_LIBRARY_DIR}/libboost_program_options.so.${_acts_boost_recommended_version} DESTINATION ${CMAKE_INSTALL_LIBDIR} OPTIONAL)

# Target for Boost::unit_test_framework

ExternalProject_Add(
  boost_build_unit_test_framework
  DOWNLOAD_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ./b2 ${Boost_B2_FLAGS} --with-test
  INSTALL_COMMAND ""
  DEPENDS boost_bootstrap
  SOURCE_DIR ${Boost_SOURCE_DIR}
  LOG_BUILD true
  BUILD_IN_SOURCE true
  BUILD_BYPRODUCTS ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.so.${_acts_boost_recommended_version}
  EXCLUDE_FROM_ALL true
)

add_library(boost_unit_test_framework SHARED IMPORTED GLOBAL)
target_link_libraries(boost_unit_test_framework INTERFACE Boost::boost)
set_target_properties(boost_unit_test_framework PROPERTIES IMPORTED_LOCATION ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.so.${_acts_boost_recommended_version})
add_dependencies(boost_unit_test_framework boost_build_unit_test_framework)
add_library(Boost::unit_test_framework ALIAS boost_unit_test_framework)
install(FILES ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.so.${_acts_boost_recommended_version} DESTINATION ${CMAKE_INSTALL_LIBDIR} OPTIONAL)
