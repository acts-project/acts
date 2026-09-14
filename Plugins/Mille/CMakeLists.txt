acts_add_library(
    PluginMille
    src/Helpers.cpp
    src/ActsToMille.cpp
    src/MillePedeSolver.cpp
    src/MillePedeError.cpp
    src/runChildProcess.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginMille
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# The solver runs `pede` as a child process via Boost.Process: the compiled
# library from Boost 1.88 on, the header-only v1 API (which builds on
# Boost.Filesystem) before that. Both are probed as optional components at the
# top level, since only this plugin needs them.
if(Boost_VERSION VERSION_GREATER_EQUAL 1.88.0)
    set(_mille_boost_component process)
    set(_mille_boost_reason "which is a compiled library since Boost 1.88")
else()
    set(_mille_boost_component filesystem)
    set(_mille_boost_reason "whose pre-1.88 API is built on Boost.Filesystem")
endif()

if(NOT TARGET Boost::${_mille_boost_component})
    foreach(_component process filesystem)
        if(TARGET Boost::${_component})
            string(APPEND _mille_boost_found " ${_component}")
        endif()
    endforeach()
    if(NOT _mille_boost_found)
        set(_mille_boost_found " none of them")
    endif()

    message(
        FATAL_ERROR
        "The Mille plugin runs `pede` as a child process via Boost.Process, "
        "${_mille_boost_reason}, so it needs the Boost "
        "${_mille_boost_component} library. The Boost ${Boost_VERSION} "
        "installation in ${Boost_DIR} provides:${_mille_boost_found}. Please "
        "enable ${_mille_boost_component} in the Boost build (with spack: "
        "`boost +${_mille_boost_component}`), or set "
        "ACTS_BUILD_PLUGIN_MILLE=OFF."
    )
endif()

# Boost is only used by the solver implementation, not by any public header
target_link_libraries(
    ActsPluginMille
    PUBLIC ActsAlignment Mille2
    PRIVATE Boost::boost Boost::${_mille_boost_component}
)

acts_compile_headers(PluginMille GLOB include/**/*.hpp)
