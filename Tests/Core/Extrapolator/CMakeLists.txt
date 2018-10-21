add_executable (NavigatorTests NavigatorTests.cpp)
target_link_libraries (NavigatorTests PRIVATE ActsCore ActsTestsCommon)
add_test(NAME NavigatorUnitTests COMMAND NavigatorTests)

add_executable (MaterialCollectionTests MaterialCollectionTests.cpp)
target_link_libraries (MaterialCollectionTests PRIVATE ActsCore ActsTestsCommon)
add_test(NAME MaterialCollectionUnitTests COMMAND MaterialCollectionTests)

add_executable (SurfaceProviderTests SurfaceProviderTests.cpp)
target_link_libraries (SurfaceProviderTests PRIVATE ActsCore)
add_test(NAME SurfaceProviderUnitTests COMMAND SurfaceProviderTests)
  
add_executable (ExtrapolatorTests ExtrapolatorTests.cpp)
target_link_libraries (ExtrapolatorTests PRIVATE ActsCore ActsTestsCommon)
add_test(NAME ExtrapolatorUnitTests COMMAND ExtrapolatorTests)



