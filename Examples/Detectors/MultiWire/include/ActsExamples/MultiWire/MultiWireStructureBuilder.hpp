// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

/*
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"*/
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Detector/interface/IExternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/MultiWire/MultiWireHelper.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace ActsExamples {

class MultiWireStructureBuilder {
 public:
  /// @brief Configuration struct for the MultiWireStructure Builder

  struct Config {
    // The name of the detector volume component
    std::string name = "";

    // The surfaces of the Multi Wire
    std::vector<std::shared_ptr<Acts::Surface>> strawSurfaces;

    // The names of the sensitive surfaces
    std::vector<std::string> sensitiveNames;

    // The names of the sensitive surfaces
    std::vector<std::string> passiveNames;

    // The binning of the layer structure
    std::vector<Acts::Experimental::LayerStructureBuilder::Binning> lbinning;

    // a tolerance config
    float toleranceOverlap = 10.;
  };

  /// Constructor
  /// @param config The configure of the MultiWireStructureBuilder

  MultiWireStructureBuilder(const Config& config);

  ~MultiWireStructureBuilder() = default;

  /// Construct the detector component

  /// @param roots The root volumes
  /// @param gctx The Geometry Context of the current geometry
  /// @return a detector component object with the detector volume of the multilayer

  Acts::Experimental::DetectorComponent construct(
      Acts::Experimental::RootDetectorVolumes& roots,
      const Acts::GeometryContext& gctx);

 private:
  Config mCfg;

  std::array<std::pair<float, float>, 3> multiWireBounds;
};
}  // namespace ActsExamples

// namespace Acts