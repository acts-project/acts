if(ACTS_BUILD_DOCS)
    # Some of the sympy code generators are written as jupytext notebooks: the
    # build runs them as scripts, and a notebook viewer opens the very same file
    # as an explained derivation. Execute those and render them into markdown
    # pages so the explanation, and the outputs it produces, are published with
    # the rest of the documentation. See docs/render_notebook.py.
    include(ActsEnsureUv)

    set(DOXYGEN_NOTEBOOK_DIR ${CMAKE_CURRENT_BINARY_DIR}/notebooks)
    # Created unconditionally: doxygen warns about an INPUT path that does not
    # exist, and warnings are errors in CI.
    file(MAKE_DIRECTORY ${DOXYGEN_NOTEBOOK_DIR})

    # Add one notebook to the documentation. PAGE_ID is both the name of the
    # generated markdown file and the doxygen page label, so other pages can
    # refer to it with @ref.
    function(acts_render_notebook)
        set(oneValueArgs SOURCE PAGE_ID)
        cmake_parse_arguments(PARSE_ARGV 0 ARGS "" "${oneValueArgs}" "")

        set(_output ${DOXYGEN_NOTEBOOK_DIR}/${ARGS_PAGE_ID}.md)

        # The renderer brings jupytext and nbconvert itself (inline script
        # metadata); what has to be added here is whatever the notebook imports.
        # The python version tracks codegen/manifest.json, so that a notebook
        # renders under the same interpreter the build generates code with.
        add_custom_command(
            OUTPUT ${_output}
            COMMAND
                ${uv_exe} run --quiet --python 3.13 --no-project
                --with-requirements
                ${PROJECT_SOURCE_DIR}/codegen/requirements.txt --with
                ${PROJECT_SOURCE_DIR}/codegen
                ${CMAKE_CURRENT_SOURCE_DIR}/render_notebook.py ${ARGS_SOURCE}
                ${_output} --page-id ${ARGS_PAGE_ID}
            DEPENDS
                ${ARGS_SOURCE}
                ${CMAKE_CURRENT_SOURCE_DIR}/render_notebook.py
                ${PROJECT_SOURCE_DIR}/codegen/src/codegen/sympy_common.py
            COMMENT "Executing and rendering notebook ${ARGS_PAGE_ID}"
            VERBATIM
        )

        set_property(GLOBAL APPEND PROPERTY ACTS_DOCS_NOTEBOOKS ${_output})
    endfunction()

    acts_render_notebook(
        SOURCE
            ${PROJECT_SOURCE_DIR}/Core/src/Propagator/detail/generate_sympy_jac.py
        PAGE_ID codegen_sympy_jac
    )

    get_property(_docs_notebooks GLOBAL PROPERTY ACTS_DOCS_NOTEBOOKS)
    add_custom_target(docs-notebooks DEPENDS ${_docs_notebooks})

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

    add_dependencies(docs docs-notebooks)

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
        examples/material_designation.cpp
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
