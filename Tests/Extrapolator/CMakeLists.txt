add_executable (NavigatorTests NavigatorTests.cpp)
target_link_libraries (NavigatorTests PRIVATE Acts::Core)
add_test(NAME NavigatorUnitTests COMMAND NavigatorTests)

add_executable (MaterialCollectionTests MaterialCollectionTests.cpp)
target_link_libraries (MaterialCollectionTests PRIVATE Acts::Core)
add_test(NAME MaterialCollectionUnitTests COMMAND MaterialCollectionTests)
  
add_executable (ExtrapolatorTests ExtrapolatorTests.cpp)
target_link_libraries (ExtrapolatorTests PRIVATE Acts::Core)
add_test(NAME ExtrapolatorUnitTests COMMAND ExtrapolatorTests)

