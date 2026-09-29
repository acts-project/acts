// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

namespace traccc::device {

/// @brief Identity projector functor
///
/// This comes in handy when using an SoA collection with some code that was
/// originally designed for AoS collections.
///
struct identity_projector {
  /// Return the input value unchanged
  template <typename T>
  TRACCC_HOST_DEVICE constexpr T operator()(const T& value) const noexcept {
    return value;
  }

};  // struct identity_projector

}  // namespace traccc::device
