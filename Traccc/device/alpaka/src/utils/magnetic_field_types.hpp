/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/bfield/magnetic_field_types.hpp"
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
#include <covfie/cuda/backend/primitive/cuda_device_array.hpp>
#elif defined(ALPAKA_ACC_GPU_HIP_ENABLED)
#include <covfie/hip/backend/primitive/hip_device_array.hpp>
#elif defined(ALPAKA_ACC_SYCL_ENABLED)
#include <covfie/sycl/backend/primitive/sycl_device_array.hpp>
#endif

namespace traccc::alpaka {

/// Inhomogeneous B-field backend type for Alpaka

#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
/// Inhomogeneous B-field backend type using CUDA global memory
template <typename scalar_t>
using inhom_bfield_backend_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::clamp<
        covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                                 covfie::backend::cuda_device_array<
                                     covfie::vector::vector_d<scalar_t, 3>>>>>>;
#elif defined(ALPAKA_ACC_GPU_HIP_ENABLED)
template <typename scalar_t>
using inhom_bfield_backend_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::clamp<
        covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                                 covfie::backend::hip_device_array<
                                     covfie::vector::vector_d<scalar_t, 3>>>>>>;
#elif defined(ALPAKA_ACC_SYCL_ENABLED)
template <typename scalar_t>
using inhom_bfield_backend_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::clamp<
        covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                                 covfie::backend::sycl_device_array<
                                     covfie::vector::vector_d<scalar_t, 3>>>>>>;
#else
template <typename scalar_t>
using inhom_bfield_backend_t = traccc::host::inhom_bfield_backend_t<scalar_t>;
#endif

// Test that the type is a valid backend for a field
static_assert(
    covfie::concepts::field_backend<inhom_bfield_backend_t<float>>,
    "alpaka::inhom_bfield_backend_t is not a valid field backend type");
/// @brief the standard list of Alpaka bfield types to support
template <typename scalar_t>
using bfield_type_list = std::tuple<const_bfield_backend_t<scalar_t>,
                                    host::inhom_bfield_backend_t<scalar_t>,
                                    inhom_bfield_backend_t<scalar_t>>;

}  // namespace traccc::alpaka
