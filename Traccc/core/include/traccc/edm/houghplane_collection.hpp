// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/edm/container.hpp>

// System include(s).
#include <array>
#include <compare>
#include <limits>

namespace traccc::edm {

/// Interface for the @c traccc::edm::houghplane_collection class.
template <typename BASE>
class houghplane : public BASE {
 public:
  /// @name Functions inherited from the base class
  /// @{

  /// Inherit the base class's constructor(s)
  using BASE::BASE;
  /// Inherit the base class's assignment operator(s).
  using BASE::operator=;

  /// @}

  /// @name Hough Plane Information
  /// @{

  TRACCC_HOST_DEVICE
  auto& plane() { return BASE::template get<0>(); }
  
  TRACCC_HOST_DEVICE
  const auto& plane() const { return BASE::template get<0>(); }

  TRACCC_HOST_DEVICE
  auto& edges_x() { return BASE::template get<1>(); }

  TRACCC_HOST_DEVICE
  const auto& edges_x() const { return BASE::template get<1>(); }

  TRACCC_HOST_DEVICE
  auto& edges_y() { return BASE::template get<2>(); }
  
  TRACCC_HOST_DEVICE
  const auto& edges_y() const { return BASE::template get<2>(); }

};
 
/// SoA container describing reconstructed hough planes
template<unsigned int N_COLS, 
	 unsigned int N_ROWS>
using houghplane_collection =
    vecmem::edm::container<houghplane,
                           // bins of the plane
                           vecmem::edm::type::vector<std::array<unsigned long long, N_COLS * N_ROWS>>,
                           // edges (x) of the planes
                           vecmem::edm::type::scalar<std::array<unsigned int, N_COLS>>,
                           // edges (y) of the planes
                           vecmem::edm::type::scalar<std::array<unsigned int, N_ROWS>>>;

}  // namespace traccc::edm

