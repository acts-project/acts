acts_add_library(
    PluginFastJet
    SHARED
    src/TrackJets.cpp
    ACTS_INCLUDE_FOLDER include/Acts
)

target_include_directories(
    ActsPluginFastJet
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginFastJet PUBLIC Acts::Core FastJet::FastJet)
