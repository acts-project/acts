# This macro adds the plugin binding as a separate library if configured to be built.
# The plugin is then added to the _plugins_built list which is propagated to the parent scope
# and used to populate the __init__.py file automatically.
macro(add_plugin_binding_if _name _dependency _flag)
    string(TOLOWER ${_name} _name_lower)
    set(_python_install_dir "${ACTS_PYTHON_INSTALL_PREFIX}/acts")

    if(${_flag})
        # automatically populate the init file with actually built modules
        list(APPEND _plugins_built "'${_name_lower}'")
        # automatically prefix the target name
        set(_target "ActsPluginsPythonBindings${_name}")
        pybind11_add_module(${_target} src/${_name}.cpp)

        install(TARGETS ${_target} DESTINATION ${_python_install_dir})

        set_target_properties(
            ${_target}
            PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
        )
        set_target_properties(
            ${_target}
            PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
        )

        target_link_libraries(
            ${_target}
            PUBLIC Acts::Core Acts::PythonUtilities ${_dependency}
        )
        add_dependencies(ActsPythonBindings ${_target})
    endif()

    set(PLUGIN_NAME "${_name}")
    set(PLUGIN_BUILD_FLAG "${_flag}")

    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/python/plugin.py.in
        ${_python_dir}/acts/${_name_lower}.py
        @ONLY
    )
    install(
        FILES ${_python_dir}/acts/${_name_lower}.py
        DESTINATION ${_python_install_dir}
    )
endmacro()

add_plugin_binding_if(Covfie Acts::PluginCovfie ACTS_BUILD_PLUGIN_TRACCC)
add_plugin_binding_if(Detray Acts::PluginDetray ACTS_BUILD_PLUGIN_TRACCC)
add_plugin_binding_if(DD4hep ActsPluginDD4hep ACTS_BUILD_PLUGIN_DD4HEP)
add_plugin_binding_if(Geant4 ActsPluginGeant4 ACTS_BUILD_PLUGIN_GEANT4)
add_plugin_binding_if(GeoModel Acts::PluginGeoModel ACTS_BUILD_PLUGIN_GEOMODEL)
add_plugin_binding_if(Hashing Acts::PluginHashing ACTS_BUILD_PLUGIN_HASHING)
add_plugin_binding_if(Json Acts::PluginJson ACTS_BUILD_PLUGIN_JSON)
add_plugin_binding_if(Svg Acts::PluginActSVG ACTS_BUILD_PLUGIN_ACTSVG)
add_plugin_binding_if(TGeo Acts::PluginRoot ACTS_BUILD_PLUGIN_ROOT)
add_plugin_binding_if(Root Acts::PluginRoot ACTS_BUILD_PLUGIN_ROOT)

# Propagate the list of built plugins to the parent scope
set(_plugins_built "${_plugins_built}" PARENT_SCOPE)
