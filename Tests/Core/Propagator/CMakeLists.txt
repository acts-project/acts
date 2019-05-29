add_executable (AbortListTests AbortListTests.cpp)
target_link_libraries (AbortListTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME AbortListUnitTest COMMAND AbortListTests)

add_executable (ActionListTests ActionListTests.cpp)
target_link_libraries (ActionListTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME ActionListUnitTests COMMAND ActionListTests)

add_executable (IntersectionCorrectorTests IntersectionCorrectorTests.cpp)
target_link_libraries (IntersectionCorrectorTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME IntersectionCorrectorUnitTests COMMAND IntersectionCorrectorTests)

add_executable (LoopProtectionTests LoopProtectionTests.cpp)
target_link_libraries (LoopProtectionTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME LoopProtectionUnitTests COMMAND LoopProtectionTests)

add_executable (ConstrainedStepTests ConstrainedStepTests.cpp)
target_link_libraries (ConstrainedStepTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME ConstrainedStepUnitTest COMMAND ConstrainedStepTests)

add_executable (EigenStepperBenchmark EigenStepperBenchmark.cpp)
target_include_directories (EigenStepperBenchmark PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (EigenStepperBenchmark PRIVATE ActsCore)
target_link_libraries (EigenStepperBenchmark PRIVATE ${Boost_LIBRARIES})

install (TARGETS EigenStepperBenchmark RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/Profiling)

add_executable (AuctioneerTests AuctioneerTests.cpp)
target_link_libraries (AuctioneerTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME AuctioneerUnitTest COMMAND AuctioneerTests)

add_executable (StepperTests StepperTests.cpp)
target_link_libraries (StepperTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME StepperUnitTest COMMAND StepperTests)

add_executable (NavigatorTests NavigatorTests.cpp)
target_link_libraries (NavigatorTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME NavigatorUnitTests COMMAND NavigatorTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST NavigatorUnitTests TARGETS NavigatorTests)

add_executable (MaterialCollectionTests MaterialCollectionTests.cpp)
target_link_libraries (MaterialCollectionTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME MaterialCollectionUnitTests COMMAND MaterialCollectionTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST MaterialCollectionUnitTests TARGETS MaterialCollectionTests)
    
add_executable (ExtrapolatorTests ExtrapolatorTests.cpp)
target_link_libraries (ExtrapolatorTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME ExtrapolatorUnitTests COMMAND ExtrapolatorTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST ExtrapolatorUnitTests TARGETS ExtrapolatorTests)

add_executable (KalmanExtrapolatorTests KalmanExtrapolatorTests.cpp)
target_link_libraries (KalmanExtrapolatorTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME KalmanExtrapolatorUnitTests COMMAND KalmanExtrapolatorTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST KalmanExtrapolatorUnitTests TARGETS KalmanExtrapolatorTests)

# these tests rely on the Legacy code
if (ACTS_BUILD_LEGACY)
  add_executable (JacobianTests JacobianTests.cpp)
  target_link_libraries (JacobianTests PRIVATE ActsCore ActsLegacy ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  add_test(NAME JacobianUnitTests COMMAND JacobianTests)
  
  add_executable (PropagatorTests PropagatorTests.cpp)
  target_link_libraries (PropagatorTests PRIVATE ActsCore ActsLegacy ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  add_test(NAME PropagatorUnitTests COMMAND PropagatorTests)
  
  # define executables
  add_executable (AtlasStepperBenchmark AtlasStepperBenchmark.cpp)
  target_include_directories (AtlasStepperBenchmark PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries (AtlasStepperBenchmark PRIVATE ActsCore ActsLegacy)
  target_link_libraries (AtlasStepperBenchmark PRIVATE ${Boost_LIBRARIES})
  
  install (TARGETS AtlasStepperBenchmark RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/Profiling)
endif()
