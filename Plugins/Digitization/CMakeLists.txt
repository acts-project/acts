file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ActsDigitizationPlugin SHARED ${src_files})

target_include_directories(ActsDigitizationPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ActsDigitizationPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ActsDigitizationPlugin PUBLIC ActsCore)

install (TARGETS ActsDigitizationPlugin EXPORT ActsDigitizationPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT ActsDigitizationPlugin TARGETS ActsDigitizationPlugin)
