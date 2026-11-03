// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_state_collection.hpp"

namespace traccc::edm {

/// Create a track state with default values.
///
/// @param measurements The collection of measurements to use for initialization
/// @param mindex       The index of the measurement to associate with the state
///
/// @return A track state object initialized with default values
///
template <typename algebra_t>
TRACCC_HOST_DEVICE
    typename track_state_collection<algebra_t>::device::object_type
    make_track_state(const measurement_collection::const_device& measurements,
                     unsigned int mindex);

}  // namespace traccc::edm

// Include the implementation.
#include "traccc/edm/impl/track_state_helpers.ipp"
