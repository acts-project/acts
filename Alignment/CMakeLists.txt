acts_add_library(
    Alignment
    SHARED
    src/Kernel/detail/AlignmentEngine.cpp
    ACTS_INCLUDE_FOLDER include/ActsAlignment
)

target_include_directories(
    ActsAlignment
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsAlignment PUBLIC Acts::Core)
install(DIRECTORY include/ActsAlignment DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
