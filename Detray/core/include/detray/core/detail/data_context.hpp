// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"

namespace detray {

namespace detail {

/// @brief ACTS-style context object.
///
/// @note the detray version simply wraps the index into the corresponding data
/// store collection.
class data_context {
 public:
  /// Default constructor
  constexpr data_context() noexcept = default;

  /// Construct from @param value .
  DETRAY_HOST_DEVICE
  explicit data_context(dindex value) : m_data{value} {}

  /// @returns the index to find the data for the context - const
  DETRAY_HOST_DEVICE
  dindex get() const { return m_data; }

 private:
  dindex m_data{0};
};

}  // namespace detail

/// Placeholder context type
class empty_context {};

/// Context type for geometry data
class geometry_context : public detail::data_context {
  using base_t = detail::data_context;

 public:
  using base_t::base_t;
};

/// Context type for magnetic field data
struct magnetic_field_context : public detail::data_context {
  using base_t = detail::data_context;

 public:
  using base_t::base_t;
};

}  // namespace detray
