// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/csv/make_measurement_hit_id_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<measurement_hit_id> make_measurement_hit_id_reader(
    std::string_view filename) {
  return {filename.data(), {"measurement_id", "hit_id"}};
}

}  // namespace traccc::io::csv
