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
#include <string>

namespace traccc::io::csv {

/// Type used in reading CSV measurement data into memory
struct measurement {
  std::uint64_t measurement_id = 0;
  std::uint64_t geometry_id = 0;
  std::uint8_t local_key = 0;
  float local0 = 0.f;
  float local1 = 0.f;
  float phi = 0.f;
  float theta = 0.f;
  float time = 0.f;
  float var_local0 = 0.f;
  float var_local1 = 0.f;
  float var_phi = 0.f;
  float var_theta = 0.f;
  float var_time = 0.f;

  DFE_NAMEDTUPLE(measurement, measurement_id, geometry_id, local_key, local0,
                 local1, phi, theta, time, var_local0, var_local1, var_phi,
                 var_theta, var_time);
};

}  // namespace traccc::io::csv
