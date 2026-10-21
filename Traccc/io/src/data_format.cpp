// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/data_format.hpp"

// System include(s).
#include <iostream>

namespace traccc {

std::ostream& operator<<(std::ostream& out, data_format format) {
  switch (format) {
    case data_format::csv:
      out << "csv";
      break;
    case data_format::binary:
      out << "binary";
      break;
    case data_format::json:
      out << "json";
      break;
    case data_format::obj:
      out << "wavefront obj";
      break;
    default:
      out << "?!?unknown?!?";
      break;
  }
  return out;
}

}  // namespace traccc
