set(unittest_extra_libraries ActsExamplesFramework)
add_unittest(Measurement MeasurementTests.cpp)
add_unittest(MuonSpacePointId MuonSpacePointIdTests.cpp)
add_unittest(JetsTests JetsTests.cpp)

#[[
if(ACTS_ENABLE_CUDA)
  add_unittest(CudaMuonPatternContainers CudaMuonPatternContainersTests.cpp)

  find_package(CUDAToolkit REQUIRED)

  target_link_libraries(
    ActsUnitTestCudaMuonPatternContainers
    PRIVATE
      CUDA::cudart
      ROOT::TreePlayer
      ROOT::Tree
      ROOT::RIO
      ROOT::Core
  )

  target_compile_definitions(
    ActsUnitTestCudaMuonPatternContainers
    # Change to public later when it becomes needed there
    PRIVATE
      ACTS_ENABLE_CUDA
  )
endif()
]]
