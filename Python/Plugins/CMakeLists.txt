set(_python_dir "${CMAKE_BINARY_DIR}/python")
set(_python_install_dir "python/acts")


pybind11_add_module(ActsPluginsPythonBindings src/PluginsModuleEntry.cpp src/Obj.cpp)

install(TARGETS ActsPluginsPythonBindings DESTINATION ${_python_install_dir}/plugins)

set_target_properties(
    ActsPluginsPythonBindings
    PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
)
set_target_properties(
    ActsPluginsPythonBindings
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
)

target_link_libraries(ActsPluginsPythonBindings PUBLIC ActsCore ActsPythonUtilities)

add_dependencies(ActsPythonBindings ActsPluginsPythonBindings)

set(py_plugin_files
    __init__.py
)

if(ACTS_BUILD_PLUGIN_ROOT)
    target_link_libraries(ActsPluginsPythonBindings PUBLIC ActsPluginRoot)
    target_sources(ActsPluginsPythonBindings PRIVATE src/TGeo.cpp)
else()
    target_sources(ActsPluginsPythonBindings PRIVATE stubs/TGeoStub.cpp)
endif()

if(ACTS_BUILD_PLUGIN_JSON)
    target_link_libraries(ActsPluginsPythonBindings PUBLIC ActsPluginJson)
    target_sources(ActsPluginsPythonBindings PRIVATE src/Json.cpp)
else()
    target_sources(ActsPluginsPythonBindings PRIVATE stubs/JsonStub.cpp)
endif()

if(ACTS_BUILD_PLUGIN_ACTSVG)
    target_link_libraries(ActsPluginsPythonBindings PUBLIC ActsExamplesIoSvg)
    target_sources(ActsPluginsPythonBindings PRIVATE src/Svg.cpp)
else()
    target_sources(ActsPluginsPythonBindings PRIVATE stubs/SvgStub.cpp)
endif()

if(ACTS_BUILD_PLUGIN_FPEMON)
    target_link_libraries(ActsPluginsPythonBindings PUBLIC ActsPluginFpeMonitoring)
    target_sources(ActsPluginsPythonBindings PRIVATE src/FpeMonitoring.cpp)
else()
    target_sources(ActsPluginsPythonBindings PRIVATE stubs/FpeMonitoringStub.cpp)
endif()

if(ACTS_BUILD_PLUGIN_GEOMODEL)
    target_link_libraries(
        ActsPluginsPythonBindings PUBLIC ActsPluginGeoModel
    )
    target_sources(ActsPluginsPythonBindings PRIVATE src/GeoModel.cpp)
else()
    target_sources(ActsPluginsPythonBindings PRIVATE stubs/GeoModelStub.cpp)
endif()

if(ACTS_BUILD_PLUGIN_TRACCC)
    target_link_libraries(ActsPluginsPythonBindings PUBLIC ActsPluginDetray)
    target_sources(ActsPluginsPythonBindings PRIVATE src/Detray.cpp)
    target_link_libraries(ActsPluginsPythonBindings PUBLIC ActsPluginCovfie)
    target_sources(ActsPluginsPythonBindings PRIVATE src/Covfie.cpp)
else()
    target_sources(ActsPluginsPythonBindings PRIVATE stubs/DetrayStub.cpp)
    target_sources(ActsPluginsPythonBindings PRIVATE stubs/CovfieStub.cpp)
    target_sources(ActsPluginsPythonBindings PRIVATE stubs/TracccStub.cpp)
endif()

# Dedicated libraries
#
# DD4hep - needs special treatment to avoid impport error
if(ACTS_BUILD_PLUGIN_DD4HEP)
    pybind11_add_module(ActsPluginsPythonBindingsDD4hep src/DD4hepComponent.cpp)
    target_link_libraries(
        ActsPluginsPythonBindingsDD4hep
        PUBLIC ActsCore ActsPythonUtilities ActsPluginDD4hep
    )
    add_dependencies(ActsPluginsPythonBindings ActsPluginsPythonBindingsDD4hep)

    install(TARGETS ActsPluginsPythonBindingsDD4hep DESTINATION ${_python_install_dir}/plugins)
    set_target_properties(
        ActsPluginsPythonBindingsDD4hep
        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    )
    set_target_properties(
        ActsPluginsPythonBindingsDD4hep
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    )
endif()

# Geant4
if(ACTS_BUILD_PLUGIN_GEANT4)
    pybind11_add_module(ActsPluginsPythonBindingsGeant4 src/Geant4Component.cpp)
    target_link_libraries(
        ActsPluginsPythonBindingsGeant4
        PUBLIC ActsCore ActsPythonUtilities ActsPluginGeant4
    )
    add_dependencies(ActsPluginsPythonBindings ActsPluginsPythonBindingsGeant4)

    install(TARGETS ActsPluginsPythonBindingsGeant4 DESTINATION ${_python_install_dir}/plugins)
    set_target_properties(
        ActsPluginsPythonBindingsGeant4
        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    )
    set_target_properties(
        ActsPluginsPythonBindingsGeant4
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    )
endif()


foreach(f ${py_plugin_files})
    set(_target ${_python_dir}/acts/plugins/${f})

    get_filename_component(_dir ${_target} DIRECTORY)
    file(MAKE_DIRECTORY ${_dir})

    file(
        CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python/${f} ${_target}
        SYMBOLIC
    )

    get_filename_component(_rel ${f} DIRECTORY)

    install(
        FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/${f}
        DESTINATION ${_python_install_dir}/plugins/${_rel}
    )
endforeach()

