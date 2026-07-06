add_unittest(EstimateTrackParamsFromSeed EstimateTrackParamsFromSeedTest.cpp)
add_unittest(HoughTransformTest HoughTransformTest.cpp)
add_unittest(AdaptiveHoughTransformTest HoughAccumulatorSectionTest.cpp)
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

    set_tests_properties(StrawLineFitterTest PROPERTIES LABELS "slow")
endif()


if(ACTS_ENABLE_CUDA)
  add_unittest(CudaHoughTransformUtils CudaHoughTransformUtilsTests.cpp)

  find_package(CUDAToolkit REQUIRED)

  target_link_libraries(
    ActsUnitTestCudaHoughTransformUtils
    PRIVATE
      ActsExamplesFramework
      ActsExamplesIoRoot
      CUDA::cudart
      ROOT::TreePlayer
      ROOT::Tree
      ROOT::RIO
      ROOT::Core
  )

  target_compile_definitions(
    ActsUnitTestCudaHoughTransformUtils
    PRIVATE
      ACTS_ENABLE_CUDA
  )
endif()
