# Remove '-Werror=null-dereference' from global flags
# It gives a weird false positive in the detray code
#
# detray/builders/grid_factory.hpp:322:41: error: null pointer dereference [-Werror=null-dereference]
#  322 |             b_values[boundary::e_min_r] = axis_spans[0UL].at(0UL);
#
#

add_library(
    ActsPluginDetray
    SHARED
    src/DetrayConversionUtils.cpp
    src/DetrayConverter.cpp
    src/DetrayGeometryConverter.cpp
    src/DetrayMaterialConverter.cpp
    src/DetraySurfaceGridsConverter.cpp
)

add_dependencies(ActsPluginDetray detray::core covfie::core vecmem::core)

target_include_directories(
    ActsPluginDetray
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(
    ActsPluginDetray
    PUBLIC
        ActsCore
        ActsPluginJson
        detray::core
        detray::core_array
        detray::io
        detray::detectors
        detray::test_utils
        vecmem::core
)

install(
    TARGETS ActsPluginDetray
    EXPORT ActsPluginDetrayTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
