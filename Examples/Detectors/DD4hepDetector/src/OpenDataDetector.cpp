// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/OpenDataDetector.hpp"

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"

namespace ActsExamples {

OpenDataDetector::OpenDataDetector(const Config& cfg,
                                   const Acts::GeometryContext& gctx)
    : DD4hepDetectorBase{cfg}, m_cfg{cfg} {
  ACTS_INFO("OpenDataDetector construct");
  construct(gctx);
}

auto OpenDataDetector::config() const -> const Config& {
  return m_cfg;
}

void OpenDataDetector::construct(const Acts::GeometryContext& gctx) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using namespace Acts::UnitLiterals;

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {0_mm, 20_mm};
  Blueprint root{cfg};

  auto volBounds = std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 1_m);
  auto vol =
      std::make_unique<TrackingVolume>(Transform3::Identity(), volBounds);

  root.addStaticVolume(std::move(vol));

  BlueprintOptions options;

  m_trackingGeometry = root.construct(options, gctx, logger());
}

}  // namespace ActsExamples
