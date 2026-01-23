acts_add_library(
    ExamplesMagneticField
    src/ScalableBFieldService.cpp
    src/ToroidalField.cpp
    src/ToroidalFieldMap.cpp
)
target_include_directories(
    ActsExamplesMagneticField
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesMagneticField
    # the ROOT libraries should be private, but if we do that then the linker
    # fails with some missing ROOT symbols.
    PUBLIC Acts::ExamplesFramework ROOT::Core ROOT::Tree
)

acts_compile_headers(ExamplesMagneticField GLOB "include/**/*.hpp")

# Add Python bindings for ToroidalField if Python bindings are enabled
if(ACTS_BUILD_EXAMPLES_PYTHON_BINDINGS)
    pybind11_add_module(acts_toroidal_field src/ToroidalFieldBindings.cpp)
    target_link_libraries(
        acts_toroidal_field
        PUBLIC ActsExamplesMagneticField Acts::Core
    )
    target_include_directories(
        acts_toroidal_field
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

    # Set output directory to make module available during build
    set_target_properties(
        acts_toroidal_field
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/acts
    )

    # Install the Python module
    install(TARGETS acts_toroidal_field DESTINATION python/acts)
endif()
