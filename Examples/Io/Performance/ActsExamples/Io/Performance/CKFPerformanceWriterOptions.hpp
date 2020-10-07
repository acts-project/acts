// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

/// Add CKFPerformanceWriter options.
///
/// @param desc The options description to add options to
void addCKFPerformanceWriterOptions(Description& desc);

/// Read CKFPerformanceWriter options to create the config.
///
/// @param variables The variables to read from
CKFPerformanceWriter::Config readCKFPerformanceWriterConfig(
    const Variables& variables);

}  // namespace Options
}  // namespace ActsExamples