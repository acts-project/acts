// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <traccc/definitions/primitives.hpp>
#include <traccc/edm/measurement_collection.hpp>
#include <traccc/edm/seed_collection.hpp>
#include <traccc/edm/spacepoint_collection.hpp>
#include <traccc/edm/track_container.hpp>
#include <traccc/edm/track_parameters.hpp>

namespace ActsExamples {

/// @name Event data exchanged between the traccc device algorithms
///
/// The device buffers live in the memory of the backend that produced them.
/// @{
using TracccHostMeasurements = traccc::edm::measurement_collection::host;
using TracccDeviceMeasurements = traccc::edm::measurement_collection::buffer;
using TracccDeviceSpacepoints = traccc::edm::spacepoint_collection::buffer;
using TracccDeviceSeeds = traccc::edm::seed_collection::buffer;
using TracccDeviceTrackParameters =
    traccc::bound_track_parameters_collection_types::buffer;
using TracccTrackContainer =
    traccc::edm::track_container<traccc::default_algebra>;
using TracccHostTracks = TracccTrackContainer::host;
using TracccDeviceTracks = TracccTrackContainer::buffer;
/// @}

}  // namespace ActsExamples
