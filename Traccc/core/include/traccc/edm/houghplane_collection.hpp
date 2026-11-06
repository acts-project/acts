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

};
 
/// SoA container describing reconstructed hough planes
template<typename T>
using houghplane_collection =
    vecmem::edm::container<houghplane,
                           // bins of the plane
                           vecmem::edm::type::jagged_vector<T>>;

}  // namespace traccc::edm

