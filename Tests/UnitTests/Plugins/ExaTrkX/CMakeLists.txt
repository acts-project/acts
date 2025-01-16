set(unittest_extra_libraries ActsPluginExaTrkX ${TORCH_LIBRARIES})

add_unittest(ExaTrkXTensorConversion ExaTrkXTensorConversionTests.cpp)
add_unittest(ExaTrkXEdgeBuilding ExaTrkXEdgeBuildingTest.cpp)
add_unittest(ExaTrkXBoostTrackBuilding ExaTrkXBoostTrackBuildingTests.cpp)
add_unittest(ExaTrkXMetricHookTests ExaTrkXMetricHookTests.cpp)
if(ACTS_EXATRKX_ENABLE_CUDA)
    add_unittest(ConnectedComponentsCuda ConnectedComponentCudaTests.cu)
    # TODO this is not nice, not sure why it does not propagate through...
    set_property(
        TARGET ActsUnitTestConnectedComponentsCuda
        PROPERTY CUDA_STANDARD 20
    )
endif()
