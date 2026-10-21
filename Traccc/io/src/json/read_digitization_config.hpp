// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/io/digitization_config.hpp"

// System include(s).
#include <string_view>

namespace traccc::io::json {

/// Read the detector digitization configuration from a JSON input file
///
/// @param filename The name of the file to read the data from
/// @return An object describing the digitization configuration of the
///         detector
///
digitization_config read_digitization_config(std::string_view filename);

}  // namespace traccc::io::json
