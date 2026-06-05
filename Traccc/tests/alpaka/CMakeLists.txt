# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2022-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

include(traccc-alpaka-functions)
traccc_enable_language_alpaka()

if(alpaka_ACC_GPU_CUDA_ENABLE)
  list(APPEND DEVICE_LIBRARIES vecmem::cuda)
elseif(alpaka_ACC_GPU_HIP_ENABLE)
  list(APPEND DEVICE_LIBRARIES vecmem::hip)
elseif(alpaka_ACC_SYCL_ENABLE)
  include( traccc-compiler-options-sycl )
  list(APPEND DEVICE_LIBRARIES vecmem::sycl)
endif()

traccc_add_alpaka_test( alpaka
   alpaka_basic.cpp
   test_cca.cpp
   LINK_LIBRARIES
   GTest::gtest_main
   traccc_tests_common
   alpaka::alpaka
   vecmem::core
   traccc::alpaka
   ${DEVICE_LIBRARIES}
)
