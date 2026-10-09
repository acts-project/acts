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
#include <covfie/core/vector.hpp>

namespace traccc {

/// Constant magnetic field backend type
template <typename scalar_t>
using const_bfield_backend_t =
    ::covfie::backend::constant<::covfie::vector::vector_d<scalar_t, 3>,
                                ::covfie::vector::vector_d<scalar_t, 3>>;
// Test that the type is a valid backend for a field
static_assert(covfie::concepts::field_backend<const_bfield_backend_t<float>>,
              "const_bfield_backend_t is not a valid field backend type");

namespace host {

/// Inhomogeneous magnetic field used for IO
template <typename scalar_t>
using inhom_io_bfield_backend_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::strided<
        covfie::vector::vector_d<std::size_t, 3>,
        covfie::backend::array<covfie::vector::vector_d<scalar_t, 3>>>>>;
// Test that the type is a valid backend for a field
static_assert(covfie::concepts::field_backend<inhom_io_bfield_backend_t<float>>,
              "inhom_io_bfield_backend_t is not a valid field backend type");

/// Inhomogeneous magnetic field backend type
template <typename scalar_t>
using inhom_bfield_backend_t = covfie::backend::affine<
    covfie::backend::linear<covfie::backend::clamp<covfie::backend::strided<
        covfie::vector::vector_d<std::size_t, 3>,
        covfie::backend::array<covfie::vector::vector_d<scalar_t, 3>>>>>>;
// Test that the type is a valid backend for a field
static_assert(covfie::concepts::field_backend<inhom_bfield_backend_t<float>>,
              "inhom_bfield_backend_t is not a valid field backend type");

/// @brief The standard list of host bfield types to support
template <typename scalar_t>
using bfield_type_list = std::tuple<const_bfield_backend_t<scalar_t>,
                                    inhom_bfield_backend_t<scalar_t>>;

}  // namespace host
}  // namespace traccc
