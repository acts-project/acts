if(ACTS_BUILD_DOCS)
    include(FetchContent)
    FetchContent_Declare(
        doxygen-awesome-css
        GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css.git
        GIT_TAG 1f3620084ff75734ed192101acf40e9dff01d848
    )
    FetchContent_MakeAvailable(doxygen-awesome-css)

    # Save the location the files were cloned into
    # This allows us to get the path to doxygen-awesome.css
    FetchContent_GetProperties(doxygen-awesome-css SOURCE_DIR AWESOME_CSS_DIR)

    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    doxygen_add_docs(
        docs
        COMMENT "Generating documentation using doxygen"
        CONFIG_FILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    )

    # docs build is off by default
    set_target_properties(docs PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()

# Additional target example build that is off by default but can be build
# without the rest of the documentation
add_library(docs-examples SHARED)
set_target_properties(docs-examples PROPERTIES EXCLUDE_FROM_ALL TRUE)
target_sources(
    docs-examples
    PRIVATE
        examples/logging.cpp
        examples/units.cpp
        examples/errors.cpp
        examples/geometry_module.cpp
)
target_link_libraries(docs-examples PRIVATE Acts::Core)

if(ACTS_BUILD_PLUGIN_DD4HEP)
    target_link_libraries(docs-examples PRIVATE Acts::PluginDD4hep)
    target_compile_definitions(docs-examples PRIVATE ACTS_BUILD_PLUGIN_DD4HEP)
endif()

# Compile the module writing templates as real geometry modules so the
# snippet examples are verified at build time.
#! [Plain Module CMake]
include(ActsGeometryModuleHelpers)
acts_add_geometry_module(
    docs-geometry-module-template
    examples/geometry_module_template.cpp
)
#! [Plain Module CMake]
add_dependencies(docs-examples docs-geometry-module-template)

if(ACTS_BUILD_PLUGIN_DD4HEP)
    #! [DD4hep Module CMake]
    include(ActsDD4hepGeometryModuleHelpers)
    acts_add_dd4hep_geometry_module(
        docs-geometry-module-template-dd4hep
        examples/geometry_module_template_dd4hep.cpp
    )
    #! [DD4hep Module CMake]
    add_dependencies(docs-examples docs-geometry-module-template-dd4hep)
endif()
