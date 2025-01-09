// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"

#include <cstddef>
#include <memory>

namespace Acts {
class IMaterialDecorator;
}

namespace ActsExamples {

class GenericDetectorElement;

class GenericDetector : public Detector {
 public:
  struct Config {
    std::size_t buildLevel = 3;
    Acts::Logging::Level logLevel = Acts::Logging::INFO;
    Acts::Logging::Level surfaceLogLevel = Acts::Logging::INFO;
    Acts::Logging::Level layerLogLevel = Acts::Logging::INFO;
    Acts::Logging::Level volumeLogLevel = Acts::Logging::INFO;
    bool buildProto = false;
    std::shared_ptr<const Acts::IMaterialDecorator> materialDecorator;
  };

  explicit GenericDetector(const Config& cfg);

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
