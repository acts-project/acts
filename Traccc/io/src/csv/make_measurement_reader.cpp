// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/csv/make_measurement_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<measurement> make_measurement_reader(
    std::string_view filename) {
  return {filename.data(),
          {"measurement_id", "geometry_id", "local_key", "local0", "local1",
           "phi", "theta", "time", "var_local0", "var_local1", "var_phi",
           "var_theta", "var_time"}};
}

}  // namespace traccc::io::csv
