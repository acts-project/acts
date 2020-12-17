// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <memory>
#include <vector>

namespace ActsExamples {
namespace Telescope {
class TelescopeDetectorElement;
}  // namespace Telescope
}  // namespace ActsExamples

struct TelescopeDetector : public ActsExamples::IBaseDetector {
  using DetectorElement = ActsExamples::Telescope::TelescopeDetectorElement;
  using DetectorElementPtr = std::shared_ptr<DetectorElement>;
  using DetectorStore = std::vector<DetectorElementPtr>;

  /// The store of the detector elements (lifetime: job)
  DetectorStore detectorStore;

  void addOptions(ActsExamples::Options::Description& desc) const override;

  std::pair<ActsExamples::IBaseDetector::TrackingGeometryPtr, ContextDecorators>
  finalize(const boost::program_options::variables_map& vm,
           std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) override;
};
