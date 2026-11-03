// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/io/csv/dfe.hpp"

// System include(s).
#include <cstdint>

namespace traccc::io::csv {

/// Type used in reading CSV measurement-to-hit ID data into memory
struct measurement_hit_id {
  std::uint64_t measurement_id = 0u;
  std::uint64_t hit_id = 0u;

  // measurement_id, hit_id
  DFE_NAMEDTUPLE(measurement_hit_id, measurement_id, hit_id);
};

}  // namespace traccc::io::csv
