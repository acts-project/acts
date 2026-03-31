set(_python_dir "${CMAKE_BINARY_DIR}/python")
set(_utilities_python_install_dir "${_acts_python_install_dir}/utilities")

set(_utilities_py_files
    __init__.py
    geometry_identifier.py
    projection2d.py
    tracking_surfaces.py
)

foreach(f ${_utilities_py_files})
    set(_target ${_python_dir}/acts/utilities/${f})
    get_filename_component(_dir ${_target} DIRECTORY)
    file(MAKE_DIRECTORY ${_dir})

    file(
        CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python/${f} ${_target}
        SYMBOLIC
    )

    install(
        FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/${f}
        DESTINATION ${_utilities_python_install_dir}
    )
endforeach()

acts_add_library(PythonUtilities src/WhiteBoardRegistry.cpp)

target_include_directories(
    ActsPythonUtilities
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsPythonUtilities PUBLIC ActsCore pybind11::pybind11)

install(
    TARGETS ActsPythonUtilities
    EXPORT ActsPythonUtilitiesTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ActsPython DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
