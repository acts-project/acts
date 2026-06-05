# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include( traccc-alpaka-functions )
include( traccc-compiler-options-cpp )

set(PUBLIC_LIBRARIES traccc::core traccc::device_common vecmem::core covfie::core)
set(PRIVATE_LIBRARIES alpaka::alpaka)

traccc_enable_language_alpaka()

if(alpaka_ACC_GPU_CUDA_ENABLE)
  list(APPEND PRIVATE_LIBRARIES CUDA::cudart vecmem::cuda traccc::cuda_utils covfie::cuda)
elseif(alpaka_ACC_GPU_HIP_ENABLE)
  find_package( HIPToolkit REQUIRED )
  list(APPEND PRIVATE_LIBRARIES HIP::hiprt vecmem::hip covfie::hip)
elseif(alpaka_ACC_SYCL_ENABLE)
  include( traccc-compiler-options-sycl )
  list(APPEND PRIVATE_LIBRARIES vecmem::sycl traccc::sycl_utils covfie::sycl)
endif()

traccc_add_alpaka_library( traccc_alpaka alpaka TYPE SHARED
  # Utility definitions.
  "include/traccc/alpaka/utils/algorithm_base.hpp"
  "src/utils/algorithm_base.cpp"
  "include/traccc/alpaka/utils/get_device_info.hpp"
  "include/traccc/alpaka/utils/make_magnetic_field.hpp"
  "src/utils/make_magnetic_field.cpp"
  "src/utils/get_device_info.cpp"
  "include/traccc/alpaka/utils/queue.hpp"
  "src/utils/queue.cpp"
  "src/utils/get_queue.hpp"
  "src/utils/get_queue.cpp"
  "include/traccc/alpaka/utils/vecmem_objects.hpp"
  "src/utils/vecmem_objects.cpp"
  "src/utils/oneDPL.hpp"
  "src/utils/parallel_algorithms.hpp"
  # Clusterization
  "include/traccc/alpaka/clusterization/clusterization_algorithm.hpp"
  "src/clusterization/clusterization_algorithm.cpp"
  "include/traccc/alpaka/clusterization/measurement_sorting_algorithm.hpp"
  "src/clusterization/measurement_sorting_algorithm.cpp"
  # Seeding code
  "include/traccc/alpaka/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
  "src/seeding/silicon_pixel_spacepoint_formation_algorithm.cpp"
  "include/traccc/alpaka/seeding/triplet_seeding_algorithm.hpp"
  "src/seeding/triplet_seeding_algorithm.cpp"
  "include/traccc/alpaka/seeding/seed_parameter_estimation_algorithm.hpp"
  "src/seeding/seed_parameter_estimation_algorithm.cpp"
  # Track finding algorithm(s).
  "include/traccc/alpaka/finding/combinatorial_kalman_filter_algorithm.hpp"
  "src/finding/combinatorial_kalman_filter_algorithm.cpp"
  "src/finding/combinatorial_kalman_filter.hpp"
  # Track fitting algorithm(s).
  "include/traccc/alpaka/fitting/kalman_fitting_algorithm.hpp"
  "src/fitting/kalman_fitting_algorithm.cpp"
  "src/fitting/kalman_fitting.hpp"
)

# Set up Thrust specifically for the traccc::alpaka library.
if(alpaka_ACC_GPU_CUDA_ENABLE)
  thrust_create_target( traccc::cuda_thrust
    HOST CPP
    DEVICE CUDA )
  target_link_libraries(traccc::cuda_thrust INTERFACE CUDA::cudart)
  list(APPEND PRIVATE_LIBRARIES traccc::cuda_thrust)
elseif(alpaka_ACC_GPU_HIP_ENABLE)
  list(APPEND PRIVATE_LIBRARIES roc::rocthrust)
else()
  thrust_create_target (traccc::alpaka_thrust
    HOST CPP
    DEVICE CPP )
  list(APPEND PRIVATE_LIBRARIES traccc::alpaka_thrust)
endif()

target_link_libraries(traccc_alpaka PUBLIC ${PUBLIC_LIBRARIES} PRIVATE ${PRIVATE_LIBRARIES})
