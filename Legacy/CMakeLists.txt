file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp" "include/*.ipp")

add_library (ActsLegacy SHARED ${src_files})

target_include_directories (ActsLegacy SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
target_include_directories (ActsLegacy PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (ActsLegacy PRIVATE ActsCore)

install (TARGETS ActsLegacy
         EXPORT ActsLegacyTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         COMPONENT ActsLegacy)
install (DIRECTORY include/Acts
         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
         COMPONENT ActsLegacy)

acts_add_targets_to_cdash_project(PROJECT ALegacy TARGETS ActsLegacy)
