// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <memory>
#include <vector>

using namespace Acts::UnitLiterals;

namespace ActsExamples {
namespace Telescope {

class TelescopeDetectorElement;
class TelescopeG4DetectorConstruction;

struct TelescopeDetector : public IBaseDetector {
  using DetectorElement = ActsExamples::Telescope::TelescopeDetectorElement;
  using DetectorElementPtr = std::shared_ptr<DetectorElement>;
  using DetectorStore = std::vector<DetectorElementPtr>;

  struct Config {
    std::vector<double> positions{{0, 30, 60, 120, 150, 180}};
    std::array<double, 2> offsets{{0, 0}};
    std::array<double, 2> bounds{{25, 100}};
    double thickness{80_um};
    int surfaceType{0};
    int binValue{2};
  };

  Config config;
  /// The store of the detector elements (lifetime: job)
  DetectorStore detectorStore;

  void addOptions(ActsExamples::Options::Description& desc) const override;

  std::pair<ActsExamples::IBaseDetector::TrackingGeometryPtr, ContextDecorators>
  finalize(const boost::program_options::variables_map& vm,
           std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) override;

  std::pair<ActsExamples::IBaseDetector::TrackingGeometryPtr, ContextDecorators>
  finalize(const Config& cfg,
           std::shared_ptr<const Acts::IMaterialDecorator> mdecorator);
};

}  // namespace Telescope
}  // namespace ActsExamples
