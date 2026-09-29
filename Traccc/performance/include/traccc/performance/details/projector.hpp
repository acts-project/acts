// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Library include(s).
#include "traccc/performance/details/is_same_object.hpp"

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

namespace traccc::details {

/// Factory creating instances of "comparator objects" for a given type
///
/// This level of abstraction is necessary to be able to construct comparator
/// objects that would have extra configuration parameters over the reference
/// object and the comparison uncertainty.
///
/// @tparam TYPE The type for which a comparator object should be generated
///
template <typename TYPE>
struct projector {
  static constexpr bool exists = false;
};

template <detray::concepts::algebra algebra_t>
struct projector<traccc::bound_track_parameters<algebra_t>> {
  static constexpr bool exists = true;

  float operator()(const traccc::bound_track_parameters<algebra_t>& i) {
    return static_cast<float>(i.phi());
  }
};

}  // namespace traccc::details
