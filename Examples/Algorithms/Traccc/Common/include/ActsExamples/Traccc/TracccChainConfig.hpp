// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Traccc/DetrayStore.hpp"

#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/spacepoint_formation_algorithm.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/finding/finding_algorithm.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

namespace ActsExamples::Traccc::Common {

struct TracccChainConfig {
  std::string inputCells;
  std::string inputMeasurements;
  std::string inputSpacePoints;
  std::string inputSeeds;
  std::string outputSpacePoints;
  std::string outputSeeds;
  std::string outputTracks;

  bool enableAmbiguityResolution;

  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
  std::shared_ptr<const DetrayHostStore> detrayStore = nullptr;
  std::shared_ptr<const Acts::ConstantBField> field = nullptr;
  Acts::GeometryHierarchyMap<DigiComponentsConfig> digitizationConfigs;

  traccc::seedfinder_config seedfinderConfig;
  traccc::spacepoint_grid_config spacepointGridConfig{seedfinderConfig};
  traccc::seedfilter_config seedfilterConfig;
  traccc::finding_config findingConfig;
  traccc::fitting_config fittingConfig;
  traccc::greedy_ambiguity_resolution_algorithm::config_t
      ambiguityResolutionConfig;
};

}  // namespace ActsExamples::Traccc::Common
