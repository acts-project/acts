// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/container.hpp"

namespace traccc {

/// Optional strip shape accompanying a spacepoint at the same output index.
/// Directions point from the first endpoint to the second (opposite to the
/// full-length direction used internally for strip intersection).
struct strip_spacepoint_shape {
  scalar top_half_length{0};
  scalar bottom_half_length{0};
  vector3 top_direction{};
  vector3 bottom_direction{};
  point3 top_center{};
  vector3 center_distance{};
};
using strip_spacepoint_shape_collection_types =
    collection_types<strip_spacepoint_shape>;

}  // namespace traccc
