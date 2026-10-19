// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/csv/make_surface_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<surface> make_surface_reader(
    std::string_view filename) {
  return {filename.data(),
          {"geometry_id", "cx", "cy", "cz", "rot_xu", "rot_xv", "rot_xw",
           "rot_zu", "rot_zv", "rot_zw"}};
}

}  // namespace traccc::io::csv
