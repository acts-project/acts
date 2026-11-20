add_subdirectory_if(Hashing ACTS_BUILD_EXAMPLES_HASHING)

add_unittest(MuonHoughSeeder MuonHoughSeederTests.cpp)
target_link_libraries(
    ActsUnitTestMuonHoughSeeder
    PRIVATE ActsExamplesTrackFinding ActsExamplesFramework
)

add_unittest(ExamplesChi2CkfTrackStateCreator Chi2CkfTrackStateCreatorTests.cpp)
target_link_libraries(
    ActsUnitTestExamplesChi2CkfTrackStateCreator
    PRIVATE ActsExamplesTrackFinding ActsExamplesFramework
)
