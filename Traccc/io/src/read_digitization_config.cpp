// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/read_digitization_config.hpp"

#include "json/read_digitization_config.hpp"
#include "traccc/io/utils.hpp"

// System include(s).
#include <stdexcept>
#include <string>

namespace traccc::io {

digitization_config read_digitization_config(std::string_view filename,
                                             data_format format) {
  // Construct the full filename.
  std::string full_filename = get_absolute_path(filename);

  // Decide how to read the file.
  switch (format) {
    case data_format::json:
      return json::read_digitization_config(full_filename);
    default:
      throw std::invalid_argument("Unsupported data format");
  }
}

}  // namespace traccc::io
