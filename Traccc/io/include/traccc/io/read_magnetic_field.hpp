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
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/utils/logging.hpp"

// System include(s).
#include <string_view>

namespace traccc::io {

/// Read a magnetic field from a file
///
/// @param[out] bfield     The field to fill
/// @param[in]  filename   The name of the file to read
/// @param[in]  format     The data format of the input file
/// @param[in]  logger     Logger to use for output messages
///
void read_magnetic_field(
    magnetic_field& bfield, std::string_view filename,
    data_format format = data_format::binary,
    std::unique_ptr<const Logger> logger = getDummyLogger().clone());

}  // namespace traccc::io
