// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Io/Json/JsonSurfacesWriter.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

// Add common Json writer options.
void addJsonWriterOptions(Description& desc);

/// Read the Json tracking geometry writer config.
ActsExamples::JsonSurfacesWriter::Config readJsonSurfacesWriterConfig(
    const Variables& vm);

}  // namespace Options
}  // namespace ActsExamples
