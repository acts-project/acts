file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

acts_add_library (DigitizationPlugin SHARED ${src_files})

target_include_directories(DigitizationPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(DigitizationPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(DigitizationPlugin PUBLIC Acts::Core)

install (TARGETS DigitizationPlugin EXPORT ActsDigitizationTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT DigitizationPlugin TARGETS DigitizationPlugin)