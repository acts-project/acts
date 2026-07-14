set(unittest_extra_libraries ActsExamplesFramework)
add_unittest(Measurement MeasurementTests.cpp)
add_unittest(MuonSpacePointId MuonSpacePointIdTests.cpp)
add_unittest(JetsTests JetsTests.cpp)

if(ACTS_ENABLE_CUDA)
    add_unittest(CudaMuonSpacePoint CudaMuonSpacePointTests.cpp)

    find_package(CUDAToolkit REQUIRED)

    target_link_libraries(ActsUnitTestCudaMuonSpacePoint PRIVATE CUDA::cudart)

    target_compile_definitions(
        ActsUnitTestCudaMuonSpacePoint
        PRIVATE ACTS_ENABLE_CUDA
    )
endif()
