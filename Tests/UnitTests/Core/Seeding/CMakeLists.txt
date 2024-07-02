add_executable(ActsUnitTestSeedFinder SeedFinderTest.cpp)
target_link_libraries(ActsUnitTestSeedFinder PRIVATE ActsCore Boost::boost Annoy::Annoy)

add_unittest(EstimateTrackParamsFromSeed EstimateTrackParamsFromSeedTest.cpp)
add_unittest(BinnedGroupTest BinnedGroupTest.cpp)
add_unittest(AnnoyTest AnnoyTest.cpp)
