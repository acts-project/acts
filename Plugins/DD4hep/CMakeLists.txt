acts_add_library(
    PluginDD4hep
    src/GeometryModuleHelper.cpp
    src/GeometryModuleLoader.cpp
    src/ConvertDD4hepDetector.cpp
    src/DD4hepBinningHelpers.cpp
    src/DD4hepMaterialHelpers.cpp
    src/DD4hepDetectorElement.cpp
    src/DD4hepDetectorSurfaceFactory.cpp
    src/DD4hepLayerBuilder.cpp
    src/DD4hepVolumeBuilder.cpp
    src/DD4hepFieldAdapter.cpp
    src/BlueprintBuilder.cpp
    src/OpenDataDetectorBuilder.cpp
    src/OpenDataDetectorBuilderTGeo.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

set(_acts_dd4hep_geometry_module_abi_tag
    "${ACTS_GEOMETRY_MODULE_ABI_TAG}|dd4hep-${DD4hep_VERSION}"
)

target_include_directories(
    ActsPluginDD4hep
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginDD4hep PUBLIC Acts::PluginRoot)
target_compile_definitions(
    ActsPluginDD4hep
    PRIVATE
        ACTS_DD4HEP_GEOMETRY_MODULE_ABI_TAG="${_acts_dd4hep_geometry_module_abi_tag}"
)

acts_compile_headers(PluginDD4hep GLOB include/**/*.hpp)

if(${DD4hep_VERSION} VERSION_LESS 1.11)
    target_include_directories(ActsPluginDD4hep PUBLIC ${DD4hep_INCLUDE_DIRS})
    target_link_libraries(ActsPluginDD4hep PUBLIC ${DD4hep_DDCORE_LIBRARY})
else()
    target_link_libraries(ActsPluginDD4hep PUBLIC DD4hep::DDCore DD4hep::DDRec)
endif()

# DDG4 actions, loaded by name into DD4hep simulations such as ddsim. Built
# whenever DD4hep comes with DDG4 and the Geant4 plugin is enabled.
if(ACTS_BUILD_PLUGIN_GEANT4 AND TARGET DD4hep::DDG4)
    message(STATUS "Build DDG4 actions in 'PluginDD4hepGeant4'")
    acts_add_library(
        PluginDD4hepGeant4
        SHARED
        src/DDG4/ActsMaterialRecorder.cpp
    )
    target_link_libraries(
        ActsPluginDD4hepGeant4
        PRIVATE
            Acts::Core
            Acts::PluginGeant4
            Acts::PluginRoot
            DD4hep::DDG4
            ROOT::Tree
    )

    # DDG4 finds actions through a .components file next to the library. This
    # is what dd4hep_generate_rootmap does, except that it assumes the library
    # is built into the current binary directory.
    if(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
        set(_ddg4_lib_dir "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
    else()
        set(_ddg4_lib_dir "${CMAKE_CURRENT_BINARY_DIR}")
    endif()
    set(_ddg4_components
        "${_ddg4_lib_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}ActsPluginDD4hepGeant4.components"
    )
    string(
        JOIN ":"
        _ddg4_ld_path
        "$<TARGET_FILE_DIR:ActsPluginDD4hepGeant4>"
        "$<TARGET_FILE_DIR:DD4hep::DD4hepGaudiPluginMgr>"
        "$ENV{LD_LIBRARY_PATH}"
    )
    add_custom_command(
        OUTPUT ${_ddg4_components}
        DEPENDS ActsPluginDD4hepGeant4
        COMMAND
            ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=${_ddg4_ld_path}"
            $<TARGET_FILE:DD4hep::listcomponents> -o ${_ddg4_components}
            $<TARGET_FILE_NAME:ActsPluginDD4hepGeant4>
        WORKING_DIRECTORY ${_ddg4_lib_dir}
        COMMENT "Generating DDG4 components file for ActsPluginDD4hepGeant4"
        VERBATIM
    )
    add_custom_target(
        ActsPluginDD4hepGeant4Components
        ALL
        DEPENDS ${_ddg4_components}
    )
    install(FILES ${_ddg4_components} DESTINATION ${CMAKE_INSTALL_LIBDIR})
else()
    message(STATUS "Skip DDG4 actions in 'PluginDD4hepGeant4'")
endif()
