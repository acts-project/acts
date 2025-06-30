add_library(ActsPluginFastJet SHARED src/TrackJets.cpp)

target_include_directories(
    ActsPluginFastJet
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginFastJet PUBLIC ActsCore FastJet)

install(
    TARGETS ActsPluginFastJet
    EXPORT ActsPluginFastJetTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
