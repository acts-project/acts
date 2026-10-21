// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/container.hpp"

// System include(s).
#include <cstdint>

namespace traccc {

/// Definition of a truth particle
struct particle {
  std::uint64_t particle_id;
  int particle_type;
  int process;
  point3 vertex;
  scalar time;
  vector3 momentum;
  scalar mass;
  scalar charge;
};

inline bool operator<(const particle& lhs, const particle& rhs) {
  if (lhs.particle_id < rhs.particle_id) {
    return true;
  }
  return false;
}

/// Declare all particle collection types
using particle_collection_types = collection_types<particle>;
/// Declare all particle container types
using particle_container_types = container_types<particle, unsigned int>;

}  // namespace traccc
