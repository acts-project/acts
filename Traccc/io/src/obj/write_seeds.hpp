// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"

// System include(s).
#include <string_view>

namespace traccc::io::obj {

/// Write a seed collection into a Wavefront OBJ file.
///
/// @param filename is the name of the output file
/// @param seeds is the seed collection to write
/// @param spacepoints is the spacepoint collection that the seeds reference
///
void write_seeds(std::string_view filename,
                 edm::seed_collection::const_view seeds,
                 edm::spacepoint_collection::const_view spacepoints);

}  // namespace traccc::io::obj
