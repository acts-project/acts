set(unittest_extra_libraries ActsExamplesIoRoot)

add_unittest(RootSimhitReaderWriter SimhitReaderWriterTests.cpp)


if(ACTS_ENABLE_CUDA)
    add_unittest(RootCudaMuonSpacePointReader RootCudaMuonSpacePointReaderTests.cpp)

  find_package(CUDAToolkit REQUIRED)

  target_link_libraries(
    ActsUnitTestRootCudaMuonSpacePointReader
    PRIVATE
      CUDA::cudart
      ROOT::TreePlayer
      ROOT::Tree
      ROOT::RIO
      ROOT::Core
  )

  target_compile_definitions(
    ActsUnitTestRootCudaMuonSpacePointReader
    # Change to public later when it becomes needed there
    PRIVATE
      ACTS_ENABLE_CUDA
  )
endif()
