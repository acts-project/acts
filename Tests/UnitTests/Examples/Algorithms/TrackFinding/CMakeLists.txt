add_subdirectory_if(Hashing ACTS_BUILD_EXAMPLES_HASHING)

add_unittest(MuonHoughSeeder MuonHoughSeederTests.cpp)
target_link_libraries(
    ActsUnitTestMuonHoughSeeder
    PRIVATE ActsExamplesTrackFinding ActsExamplesFramework
)

add_unittest(GlobalPatternFinderAlgorithm GlobalPatternFinderAlgorithmTests.cpp)
target_link_libraries(
    ActsUnitTestGlobalPatternFinderAlgorithm
    PRIVATE ActsExamplesTrackFinding ActsExamplesFramework
)
