/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Covfie include(s).
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/clamp.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/concepts.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/vector.hpp>
#include <covfie/hip/backend/primitive/hip_device_array.hpp>

#include "traccc/bfield/magnetic_field_types.hpp"

namespace traccc::hip {

/// Inhomogeneous B-field backend type using HIP global memory
template <typename scalar_t>
using inhom_global_bfield_backend_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::clamp<
        covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                                 covfie::backend::hip_device_array<
                                     covfie::vector::vector_d<scalar_t, 3>>>>>>;
// Test that the type is a valid backend for a field
static_assert(
    covfie::concepts::field_backend<inhom_global_bfield_backend_t<float>>,
    "hip::inhom_global_bfield_backend_t is not a valid field backend type");

/// @brief the standard list of HIP bfield types to support
template <typename scalar_t>
using bfield_type_list = std::tuple<const_bfield_backend_t<scalar_t>,
                                    inhom_global_bfield_backend_t<scalar_t>>;

}  // namespace traccc::hip
