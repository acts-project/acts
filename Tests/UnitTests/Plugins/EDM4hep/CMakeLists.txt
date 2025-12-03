add_unittest(ConvertTrackEDM4hep ConvertTrackEDM4hepTest.cpp)

target_link_libraries(
    ActsUnitTestConvertTrackEDM4hep
    PUBLIC Acts::PluginEDM4hep
)

set(unittest_extra_libraries Acts::PluginEDM4hep)
add_unittest(PodioTrackContainer PodioTrackContainerTest.cpp)
add_unittest(PodioTrackStateContainer PodioTrackStateContainerTest.cpp)
