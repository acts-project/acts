// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
#include <covfie/cuda/backend/primitive/cuda_device_array.hpp>
#include <covfie/cuda/backend/primitive/cuda_texture.hpp>

#include "traccc/bfield/magnetic_field_types.hpp"

namespace traccc::cuda {

/// Inhomogeneous B-field backend type using CUDA global memory
template <typename scalar_t>
using inhom_global_bfield_backend_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::clamp<
        covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                                 covfie::backend::cuda_device_array<
                                     covfie::vector::vector_d<scalar_t, 3>>>>>>;
// Test that the type is a valid backend for a field
static_assert(
    covfie::concepts::field_backend<inhom_global_bfield_backend_t<float>>,
    "cuda::inhom_global_bfield_backend_t is not a valid field backend type");

/// Inhomogeneous B-field backend type using CUDA texture memory
using inhom_texture_bfield_backend_t = covfie::backend::affine<
    covfie::backend::cuda_texture<covfie::vector::vector_d<float, 3>,
                                  covfie::vector::vector_d<float, 3>>>;
// Test that the type is a valid backend for a field
static_assert(
    covfie::concepts::field_backend<inhom_texture_bfield_backend_t>,
    "cuda::inhom_texture_bfield_backend_t is not a valid field backend type");

/// @brief the standard list of CUDA bfield types to support
template <typename scalar_t>
using bfield_type_list = std::tuple<const_bfield_backend_t<scalar_t>,
                                    inhom_global_bfield_backend_t<scalar_t>,
                                    inhom_texture_bfield_backend_t>;

}  // namespace traccc::cuda
