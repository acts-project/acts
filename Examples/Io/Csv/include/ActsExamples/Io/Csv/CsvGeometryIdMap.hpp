// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/GeometryContainers.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

/// Read the map between Athena and ACTS geometry identifiers from a CSV file
/// with the columns `athena_id` and `acts_id`.
/// @param path The CSV file
/// @return The geometry identifier map
std::shared_ptr<GeometryIdMapActsAthena> readCsvGeometryIdMap(
    const std::string& path);

}  // namespace ActsExamples
