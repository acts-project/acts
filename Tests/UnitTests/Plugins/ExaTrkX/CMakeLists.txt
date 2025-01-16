set(unittest_extra_libraries ActsPluginExaTrkX ${TORCH_LIBRARIES})

add_unittest(ExaTrkXTensorConversion ExaTrkXTensorConversionTests.cpp)
add_unittest(ExaTrkXEdgeBuilding ExaTrkXEdgeBuildingTest.cpp)
add_unittest(ExaTrkXBoostTrackBuilding ExaTrkXBoostTrackBuildingTests.cpp)
add_unittest(ExaTrkXMetricHookTests ExaTrkXMetricHookTests.cpp)
if(ACTS_EXATRKX_ENABLE_CUDA)
    add_unittest(ConnectedComponentsCuda ConnectedComponentCudaTests.cu)
    # not sure why it uses C++14 if I do not specify this...
    target_compile_features(
        ActsUnitTestConnectedComponentsCuda
        PUBLIC cxx_std_20
    )
endif()
