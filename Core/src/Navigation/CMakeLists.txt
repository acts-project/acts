target_sources(
    ActsCore
    PRIVATE
        NavigationStream.cpp
        TryAllNavigationPolicy.cpp
        SurfaceArrayNavigationPolicy.cpp
        CylinderNavigationPolicy.cpp
        MultiLayerNavigationPolicy.cpp
        MultiNavigationPolicy.cpp
)

if(ACTS_BUILD_PLUGIN_TRACCC)
    target_link_libraries(ActsCore PRIVATE detray::io)
    target_sources(ActsCore PRIVATE DetrayNavigation.cpp)
else()
    target_sources(ActsCore PRIVATE DetrayNavigationStub.cpp)
endif()
