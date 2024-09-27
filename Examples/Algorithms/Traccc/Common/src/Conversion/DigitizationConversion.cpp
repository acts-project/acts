// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/Conversion/DigitizationConversion.hpp"

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Traccc/DigitizationConfig.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"

#include <utility>
#include <vector>

namespace ActsExamples::Traccc::Common::Conversion {

Acts::TracccPlugin::DigitizationConfig tracccConfig(
    const Acts::GeometryHierarchyMap<DigiComponentsConfig>& config) {
  using ElementType = std::pair<Acts::GeometryIdentifier,
                                Acts::TracccPlugin::ModuleDigitizationConfig>;

  std::vector<ElementType> vec;

  for (std::size_t i = 0; i < config.size(); i++) {
    vec.push_back({config.idAt(i),
                   Acts::TracccPlugin::ModuleDigitizationConfig{
                       config.valueAt(i).geometricDigiConfig.segmentation}});
  }

  return Acts::TracccPlugin::DigitizationConfig(vec);
}

}  // namespace ActsExamples::Traccc::Common::Conversion
