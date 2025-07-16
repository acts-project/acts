pybind11_add_module(ActsPythonBindings
    src/ActsPythonModule.cpp
)

install(TARGETS ActsPythonBindings DESTINATION ${_python_install_dir})

set_target_properties(
    ActsPythonBindings
    PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
)
set_target_properties(
    ActsPythonBindings
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
)

target_link_libraries(
    ActsPythonBindings
    PUBLIC
        ActsCore
        ActsPythonUtilities
)

target_include_directories(
    ActsPythonBindings
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(DIRECTORY include/ActsPython DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
