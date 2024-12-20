set(unittest_extra_libraries ActsPluginExaTrkX ${TORCH_LIBRARIES})

add_unittest(ExaTrkXTensorConversion ExaTrkXTensorConversionTests.cpp)
add_unittest(ExaTrkXEdgeBuilding ExaTrkXEdgeBuildingTest.cpp)
#add_unittest(ExaTrkXBoostTrackBuilding ExaTrkXBoostTrackBuildingTests.cpp)
add_unittest(ExaTrkXMetricHookTests ExaTrkXMetricHookTests.cpp)
add_unittest(Walkthrough WalkThroughTests.cpp)
add_unittest(ModuleMap ModuleMapTests.cpp)
add_unittest(ConnectedComponentsCuda ConnectedComponentCudaTests.cu)
set_target_properties(
    ActsUnitTestConnectedComponentsCuda
    PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON
)
