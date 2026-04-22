// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/io/utils/file_handle.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/clamp.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/vector.hpp>

// System include(s)
#include <ios>
#include <iostream>
#include <stdexcept>
#include <string>

namespace detray {

namespace bfield {

/// Constant bfield (host and device)
template <typename T>
using const_bknd_t = covfie::backend::constant<covfie::vector::vector_d<T, 3>,
                                               covfie::vector::vector_d<T, 3>>;

template <typename T>
using const_field_t = covfie::field<const_bknd_t<T>>;

/// Inhomogeneous field (host)
template <typename T>
using inhom_bknd_t = covfie::backend::affine<
    covfie::backend::linear<covfie::backend::clamp<covfie::backend::strided<
        covfie::vector::vector_d<std::size_t, 3>,
        covfie::backend::array<covfie::vector::vector_d<T, 3>>>>>>;

// Test that the type is a valid backend for a field
static_assert(covfie::concepts::field_backend<inhom_bknd_t<float>>,
              "inhom_bknd_t is not a valid host field backend type");

/// Inhomogeneous field (host) for IO
template <typename T>
using inhom_bknd_io_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::strided<
        covfie::vector::vector_d<std::size_t, 3>,
        covfie::backend::array<covfie::vector::vector_d<T, 3>>>>>;

static_assert(covfie::concepts::field_backend<inhom_bknd_io_t<float>>,
              "inhom_bknd_io_t is not a valid host field backend type");

/// Inhomogeneous field with nearest neighbor (host)
template <typename T>
using inhom_bknd_nn_t =
    covfie::backend::affine<covfie::backend::nearest_neighbour<
        covfie::backend::clamp<covfie::backend::strided<
            covfie::vector::ulong3,
            covfie::backend::array<covfie::vector::vector_d<T, 3>>>>>>;

static_assert(covfie::concepts::field_backend<inhom_bknd_nn_t<float>>,
              "inhom_bknd_nn_t is not a valid host field backend type");

/// Inhomogeneous field with nearest neighbor (host) for IO
template <typename T>
using inhom_bknd_nn_io_t = covfie::backend::affine<
    covfie::backend::nearest_neighbour<covfie::backend::strided<
        covfie::vector::ulong3,
        covfie::backend::array<covfie::vector::vector_d<T, 3>>>>>;

static_assert(covfie::concepts::field_backend<inhom_bknd_nn_io_t<float>>,
              "inhom_bknd_nn_io_t is not a valid host field backend type");

template <typename T>
using inhom_field_t = covfie::field<inhom_bknd_t<T>>;

template <typename T>
using inhom_field_io_t = covfie::field<inhom_bknd_io_t<T>>;

}  // namespace bfield

/// @brief Function that reads the first 4 bytes of a potential bfield file and
/// checks that it contains data for a covfie field
inline bool check_covfie_file(const std::string& file_name) {
  // Open binary file
  io::file_handle file{file_name, std::ios_base::in | std::ios_base::binary};

  // See "covfie/lib/core/utility/binary_io.hpp"
  std::uint32_t hdr = covfie::utility::read_binary<std::uint32_t>(*file);

  // Compare to magic bytes
  return (hdr == covfie::utility::MAGIC_HEADER);
}

/// @brief function that reads a covfie field from file
template <typename bfield_t>
inline bfield_t read_bfield(const std::string& file_name) {
  if (!check_covfie_file(file_name)) {
    throw std::runtime_error("Not a valid covfie file: " + file_name);
  }

  // Open binary file
  io::file_handle file{file_name, std::ios_base::in | std::ios_base::binary};

  return bfield_t(*file);
}

/// @returns a constant covfie field constructed from the field vector @param B
template <typename T, concepts::vector3D vector3_t>
inline bfield::const_field_t<T> create_const_field(const vector3_t& B) {
  return bfield::const_field_t<T>{covfie::make_parameter_pack(
      typename bfield::const_bknd_t<T>::configuration_t{B[0], B[1], B[2]})};
}

/// @returns a constant covfie field constructed from the field vector @param B
template <typename T>
inline bfield::inhom_field_t<T> create_inhom_field() {
  return bfield::inhom_field_t<T>(read_bfield<bfield::inhom_field_io_t<T>>(
      !std::getenv("DETRAY_BFIELD_FILE") ? ""
                                         : std::getenv("DETRAY_BFIELD_FILE")));
}

}  // namespace detray
