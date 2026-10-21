// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s).
#include <iosfwd>

namespace traccc {

/// Format for an input or output file
enum data_format : int {
  csv = 0,     ///< Comma-separated values
  binary = 1,  ///< Binary format
  json = 2,    ///< JSON format
  obj = 3,     ///< Wavefront OBJ format
};

/// Printout helper for @c traccc::data_format
std::ostream& operator<<(std::ostream& out, data_format format);

}  // namespace traccc
