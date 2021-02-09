// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <vector>

namespace ActsExamples {
namespace Options {

/// Add TGeo detector options prefixed with geo-tgeo.
void addTGeoGeometryOptions(Description& desc);

/// Read the TGeo layer builder configurations from the user configuration.
std::vector<Acts::TGeoLayerBuilder::Config> readTGeoLayerBuilderConfigs(
    const Variables& vars);

}  // namespace Options
}  // namespace ActsExamples
