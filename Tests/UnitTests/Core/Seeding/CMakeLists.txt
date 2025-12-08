add_executable(ActsUnitTestSeedFinder SeedFinderTest.cpp)
target_link_libraries(ActsUnitTestSeedFinder PRIVATE Acts::Core Boost::boost)

add_unittest(EstimateTrackParamsFromSeed EstimateTrackParamsFromSeedTest.cpp)
add_unittest(BinnedGroupTest BinnedGroupTest.cpp)
#add_unittest(PathSeeder PathSeederTest.cpp)
add_unittest(HoughTransformTest HoughTransformTest.cpp)
add_unittest(UtilityFunctions UtilityFunctionsTests.cpp)
add_unittest(CandidatesForMiddleSp CandidatesForMiddleSpTests.cpp)
add_unittest(StrawLineResiduals StrawLineResidualTest.cpp)

if(ACTS_BUILD_PLUGIN_ROOT)
    add_unittest(FastStrawLineFitTests FastStrawLineFitTests.cpp)
    add_unittest(StrawLineFitterTest StrawLineFitterTest.cpp)
    add_unittest(StrawLineSeederTest StrawLineSeederTest.cpp)

    target_link_libraries(
        ActsUnitTestFastStrawLineFitTests
        PRIVATE ROOT::Tree ROOT::RIO
    )
    target_link_libraries(
        ActsUnitTestStrawLineFitterTest
        PRIVATE ROOT::Tree ROOT::RIO ROOT::Hist
    )
    target_link_libraries(
        ActsUnitTestStrawLineSeederTest
        PRIVATE ROOT::Tree ROOT::RIO
    )
endif()
