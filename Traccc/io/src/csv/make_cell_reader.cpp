// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/csv/make_cell_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<cell> make_cell_reader(std::string_view filename) {
  return {filename.data(),
          {"geometry_id", "measurement_id", "cannel0", "channel1", "timestamp",
           "value"}};
}

}  // namespace traccc::io::csv
