add_library(
    ActsPythonExamples
    INTERFACE
)

target_include_directories(
        ActsPythonExamples
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsPythonExamples INTERFACE ActsCore)

install(
    TARGETS ActsPythonExamples
    EXPORT ActsPythonExamplesTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ActsPython DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})