acts_add_library(PythonUtilities INTERFACE)

target_include_directories(
    ActsPythonUtilities
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsPythonUtilities INTERFACE ActsCore pybind11::pybind11)

install(
    TARGETS ActsPythonUtilities
    EXPORT ActsPythonUtilitiesTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ActsPython DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
