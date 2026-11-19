# This macro adds the plugin binding as a separate library if configured to be built.
# The plugin is then added to the _plugins_built list which is propagated to the parent scope
# and used to populate the __init__.py file automatically.
macro(add_plugin_binding_if _name _dependency _flag)
    string(TOLOWER ${_name} _name_lower)

    if(${_flag})
        # automatically populate the init file with actually built modules
        list(APPEND _plugins_built "'${_name_lower}'")
        # automatically prefix the target name
        set(_target "ActsPluginsPythonBindings${_name}")
        acts_add_python_module(
            ${_target}
            SOURCES src/${_name}.cpp
            OUTPUT_DIRECTORY ${_python_dir}/acts
            INSTALL_DIRECTORY ${_acts_python_install_dir}
            LIBRARY_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
        )

        target_link_libraries(
            ${_target}
            PUBLIC Acts::Core Acts::PythonUtilities ${_dependency}
        )
        add_dependencies(ActsPythonBindings ${_target})
        acts_pybind_generate_stubs(
            acts.${_name_lower}
            ${_target}
            ENABLED ${ACTS_GENERATE_PYTHON_STUBS}
            PYTHON_ROOT ${_python_dir}
            PYTHON_EXECUTABLE ${Python_EXECUTABLE}
            DEPENDS ActsPythonBindings
            IGNORE_INVALID_IDENTIFIERS "GeoModelTree::FpvConstLink"
        )
    endif()

    set(PLUGIN_NAME "${_name}")
    set(PLUGIN_BUILD_FLAG "${_flag}")

    # by default, we only expose "Config" objects as
    # kwargs-style constructors
    set(SYMBOLS_FOR_KWARGS Config)

    # plugins may add additional symbols to expose this way
    set(_extra_kwargs "${_name}_MAKE_EXTRA_KWARGS")
    if(DEFINED ${_extra_kwargs})
        list(APPEND SYMBOLS_FOR_KWARGS ${${_extra_kwargs}})
    endif()
    string(JOIN "\", \"" SYMBOLS_FOR_KWARGS_PY ${SYMBOLS_FOR_KWARGS})
    set(SYMBOLS_FOR_KWARGS_PY "[\"${SYMBOLS_FOR_KWARGS_PY}\"]")

    acts_python_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/python/plugin.py.in
        ${_python_dir}/acts/${_name_lower}.py
        ${_acts_python_install_dir}
    )
endmacro()

if(ACTS_BUILD_PLUGIN_MILLE)
    # Additional classes with kwarg-style
    # ctors for the Mille plugin
    set(Mille_MAKE_EXTRA_KWARGS
        Result
        MillePedeParameterResult
        MillePedeSteeringConfig
        MillePedeEqualityConstraint,
        MillePedeParameterSeed,
        MillePedeSolverConfig
    )
endif()

add_plugin_binding_if(Arrow Acts::PluginArrow ACTS_BUILD_PLUGIN_ARROW)
add_plugin_binding_if(Gnn Acts::PluginGnn ACTS_BUILD_PLUGIN_GNN)
add_plugin_binding_if(Covfie Acts::PluginCovfie ACTS_BUILD_PLUGIN_TRACCC)
add_plugin_binding_if(Detray Acts::PluginDetray ACTS_BUILD_PLUGIN_TRACCC)
add_plugin_binding_if(Vecmem Acts::PluginDetray ACTS_BUILD_PLUGIN_TRACCC)
add_plugin_binding_if(DD4hep ActsPluginDD4hep ACTS_BUILD_PLUGIN_DD4HEP)
add_plugin_binding_if(Geant4 ActsPluginGeant4 ACTS_BUILD_PLUGIN_GEANT4)
add_plugin_binding_if(GeoModel Acts::PluginGeoModel ACTS_BUILD_PLUGIN_GEOMODEL)
add_plugin_binding_if(Json Acts::PluginJson ACTS_BUILD_PLUGIN_JSON)
add_plugin_binding_if(Svg Acts::PluginActSVG ACTS_BUILD_PLUGIN_ACTSVG)
add_plugin_binding_if(TGeo Acts::PluginRoot ACTS_BUILD_PLUGIN_ROOT)
add_plugin_binding_if(Root Acts::PluginRoot ACTS_BUILD_PLUGIN_ROOT)
add_plugin_binding_if(Mille Acts::PluginMille ACTS_BUILD_PLUGIN_MILLE)

# Propagate the list of built plugins to the parent scope
set(_plugins_built "${_plugins_built}" PARENT_SCOPE)

if(TARGET ActsPluginsPythonBindingsDetray)
    target_link_libraries(
        ActsPluginsPythonBindingsDetray
        PRIVATE detray::python_common
    )
    add_dependencies(ActsPluginsPythonBindingsDetray DetrayPythonBindings)
    add_dependencies(ActsPluginsPythonBindingsVecmem DetrayPythonBindings)
endif()
