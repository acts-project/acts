acts_add_library(
    ExamplesAlignment
    SHARED
    src/AlignmentAlgorithm.cpp
    src/AlignmentAlgorithmFunction.cpp
)
target_include_directories(
    ActsExamplesAlignment
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesAlignment
    PUBLIC
        Acts::Core
        Acts::Alignment
        Acts::ExamplesFramework
        Acts::ExamplesMagneticField
)

install(
    TARGETS ActsExamplesAlignment
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
