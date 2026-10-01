// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/seeding/strip_spacepoint_formation_data.hpp"

namespace traccc {

/// Event-local pair of indices into the sorted measurement collection.
/// Carry the selected rule's tolerances through to the formation pass.
struct strip_pair {
  unsigned int measurement_index_1;
  unsigned int measurement_index_2;
  scalar strip_length_gap_tolerance;
  scalar strip_length_tolerance;
};

using strip_pair_collection_types = collection_types<strip_pair>;

}  // namespace traccc
