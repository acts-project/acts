file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp" "include/*.ipp")

add_library (ACTSLegacy SHARED ${src_files})

target_include_directories (ACTSLegacy SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
target_include_directories (ACTSLegacy PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (ACTSLegacy PRIVATE ACTSCore)

install (TARGETS ACTSLegacy
         EXPORT ACTSLegacyTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         COMPONENT Legacy)
install (DIRECTORY include/ACTS
         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
         COMPONENT Legacy)

acts_add_targets_to_cdash_project(PROJECT ALegacy TARGETS ACTSLegacy)
