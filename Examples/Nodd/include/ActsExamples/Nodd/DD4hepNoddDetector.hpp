// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

namespace ActsExamples {

/// Load nODD DD4hep; optionally build its complete DD4hep-backed Gen3 geometry.
class DD4hepNoddDetector final : public DD4hepDetectorBase {
 public:
  struct Config : DD4hepDetectorBase::Config {
    Config();
    bool gen3 = false;
    double navigationEnvelopeMm = 1.0;
  };

  explicit DD4hepNoddDetector(const Config& cfg);

  const Config& config() const override;

  /// JSON audit of source identities and actual converted transforms/bounds.
  std::string geometryReport() const;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
