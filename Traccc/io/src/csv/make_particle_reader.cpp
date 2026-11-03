// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/csv/make_particle_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<particle> make_particle_reader(
    std::string_view filename) {
  return {filename.data(),
          {"particle_id", "particle_type", "process", "vx", "vy", "vz", "vt",
           "px", "py", "pz", "m", "q"}};
}

}  // namespace traccc::io::csv
