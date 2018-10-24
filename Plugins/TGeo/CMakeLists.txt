file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ActsTGeoPlugin SHARED ${src_files})

if(NOT "$ENV{CMAKE_PREFIX_PATH}" STREQUAL "")
  string(REGEX MATCH "(LCG_[0-9][0-9])" _lcg_version "$ENV{CMAKE_PREFIX_PATH}")
  if (NOT "${_lcg_version}" STREQUAL "")
    string(REGEX MATCH "([A-Za-z0-9_-]*)$" _lcg_platform "$ENV{CMAKE_PREFIX_PATH}")
    if(NOT "${_lcg_platform}" STREQUAL "")
      if(${_lcg_version} STREQUAL "LCG_93" AND ${_lcg_platform} STREQUAL "x86_64-centos7-gcc7-opt")
        message(FATAL_ERROR "Compiling on ${_lcg_platform} using ${_lcg_version} with TGeoPlugin \
        enabled. There is no clean way to resolve C++ standard requirements. Bailing out")
      endif()
    endif()
  endif()
endif()



target_include_directories(ActsTGeoPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ActsTGeoPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ActsTGeoPlugin PUBLIC ActsCore)
target_link_libraries(ActsTGeoPlugin PUBLIC ROOT::Geom)
target_link_libraries(ActsTGeoPlugin PUBLIC IdentificationPlugin)

# PROBLEM: depending on LCG version, ROOT might have been compiled with 
# c++17 or not. Detect and set standard on this target
string(REPLACE " " ";" _root_cxx_flags_list ${ROOT_CXX_FLAGS})
if ("-std=c++1z" IN_LIST _root_cxx_flags_list
    OR "-std=c++17" IN_LIST _root_cxx_flags_list)

  # can our compiler do c++17?
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-std=c++17 _have_cxx17)
  if(_have_cxx17)
    # can our cmake set CXX_STANDARD to 17?
    if(${CMAKE_VERSION} VERSION_GREATER "3.8.0")
      target_compile_features(ActsTGeoPlugin PUBLIC cxx_std_17)
    else()
      message(FATAL_ERROR "ROOT compiled with C++17, but this CMake version can't set it")
    endif()
  else()
    message(SEND_ERROR "ROOT compiled with C++17, but compiler does not support that")
  endif()


endif()

install (TARGETS ActsTGeoPlugin EXPORT ActsTGeoPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT ActsTGeoPlugin TARGETS ActsTGeoPlugin)
