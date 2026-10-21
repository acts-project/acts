// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/edm/track_container.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"

// System include(s).
#include <string_view>

namespace traccc::io::obj {

/// Write a track candidate container into a Wavefront OBJ file.
///
/// @param filename is the name of the output file
/// @param tracks is the track container to write
/// @param detector is the Detray detector describing the geometry
///
void write_tracks(std::string_view filename,
                  edm::track_container<default_algebra>::const_view tracks,
                  const traccc::host_detector& detector);

}  // namespace traccc::io::obj
