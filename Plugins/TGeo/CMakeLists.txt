file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ActsTGeoPlugin SHARED ${src_files})

target_include_directories(ActsTGeoPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ActsTGeoPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ActsTGeoPlugin PUBLIC ActsCore)
target_link_libraries(ActsTGeoPlugin PUBLIC ROOT::Geom)
target_link_libraries(ActsTGeoPlugin PUBLIC IdentificationPlugin)

# PROBLEM: depending on LCG version, ROOT might have been compiled with c++17 or
# not. Detect and set standard on this target
set(_root_cxx_flags ${ROOT_CXX_FLAGS})
separate_arguments(_root_cxx_flags)
if("-std=c++1z" IN_LIST _root_cxx_flags OR "-std=c++17" IN_LIST _root_cxx_flags)
  # per-target standard reqirements are only possible in 3.8 or higher
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.8.0")
    target_compile_features(ActsTGeoPlugin PUBLIC cxx_std_17)
  else()
    message(FATAL_ERROR "ROOT compiled with C++17, but this CMake version cannot set it")
  endif()
endif()

install (TARGETS ActsTGeoPlugin EXPORT ActsTGeoPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT ActsTGeoPlugin TARGETS ActsTGeoPlugin)
