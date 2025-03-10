set(library_sources
    src/RootMaterialTrackPayload.cpp
    src/RootMaterialDecorator.cpp
)

add_library(ActsPluginRoot SHARED ${library_sources})

target_include_directories(
    ActsPluginRoot
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsPluginRoot PUBLIC ActsCore ROOT::Core  ROOT::Hist ROOT::Tree)

install(
    TARGETS ActsPluginRoot
    EXPORT ActsPluginRootTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
