// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/csv/make_hit_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<hit> make_hit_reader(std::string_view filename) {
  return {filename.data(),
          {"particle_id", "geometry_id", "tx", "ty", "tz", "tt", "tpx", "tpy",
           "tpz", "te", "deltapx", "deltapy", "deltapz", "deltae", "index"}};
}

}  // namespace traccc::io::csv
