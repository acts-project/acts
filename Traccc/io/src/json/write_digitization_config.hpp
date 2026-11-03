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

/// Write a digitization configuration to a file
///
/// @param filename The name of the file to write the data to
/// @param config The digitization configuration to write
///
void write_digitization_config(std::string_view filename,
                               const digitization_config& config);

}  // namespace traccc::io::json
