file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ActsJsonPlugin INTERFACE)

target_include_directories (ActsJsonPlugin SYSTEM INTERFACE ${Boost_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
target_include_directories (ActsJsonPlugin INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

install (TARGETS ActsJsonPlugin
         EXPORT ActsJsonPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

