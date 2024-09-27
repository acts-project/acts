add_library(ActsPluginTraccc SHARED src/CellConversion.cpp)

target_include_directories(
    ActsPluginTraccc
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsPluginTraccc
    PUBLIC
        ActsCore
        traccc::core
        vecmem::core
        detray::test_utils
        detray::io
        ActsPluginCovfie
)

install(
    TARGETS ActsPluginTraccc
    EXPORT ActsPluginTracccTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
