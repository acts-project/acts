// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/io/data_format.hpp"

// Project include(s).
#include "traccc/io/conditions_config.hpp"

// System include(s).
#include <string_view>

namespace traccc::io {

/// Read the detector digitization configuration from an input file
///
/// @param filename The name of the file to read the data from
/// @param format The format of the input file
/// @return An object describing the digitization configuration of the
///         detector
///
conditions_config read_conditions_config(
    std::string_view filename, data_format format = data_format::json);

}  // namespace traccc::io
