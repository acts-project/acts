// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

namespace Acts {
class GeometryContext;
}

namespace ActsExamples {

class OpenDataDetector final : public DD4hepDetectorBase {
 public:
  struct Config : public DD4hepDetectorBase::Config {};

  explicit OpenDataDetector(const Config& cfg,
                            const Acts::GeometryContext& gctx);

  const Config& config() const override;

 private:
  void construct(const Acts::GeometryContext& gctx);

  Config m_cfg;
};

}  // namespace ActsExamples
