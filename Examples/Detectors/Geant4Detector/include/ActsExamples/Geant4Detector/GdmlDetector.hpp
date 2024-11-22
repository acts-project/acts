// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DetectorCommons/DetectorBase.hpp"

#include <memory>

namespace ActsExamples {

struct GdmlDetectorFactory : public DetectorFactoryBase {
  struct Config {
    std::string path;
  };

  explicit GdmlDetectorFactory(const Config& cfg);

  std::shared_ptr<DetectorBase> buildDetector() const override;

 private:
  Config m_cfg;
};

class GdmlDetector : public PreConstructedDetector {
 public:
  explicit GdmlDetector(const GdmlDetectorFactory::Config& config);

  std::unique_ptr<G4VUserDetectorConstruction> buildGeant4DetectorConstruction(
      const Geant4ConstructionOptions& options) const override;

 private:
  GdmlDetectorFactory::Config m_cfg;
};

}  // namespace ActsExamples
