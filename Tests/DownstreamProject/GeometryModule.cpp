// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryModuleHelper.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

#include <memory>
#include <stdexcept>

namespace {

std::unique_ptr<Acts::TrackingGeometry> buildGeometryModule(
    const Acts::Logger& logger) {
  using namespace Acts;

  ACTS_INFO("Building geometry module");

  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  Experimental::Blueprint::Config cfg;
  cfg.envelope = ExtentEnvelope{{
      .x = {10., 10.},
      .y = {10., 10.},
      .z = {10., 10.},
  }};

  Experimental::Blueprint root{cfg};

  auto outerBounds = std::make_shared<CuboidVolumeBounds>(1000., 1000., 1000.);
  auto outerVol = std::make_unique<TrackingVolume>(Transform3::Identity(),
                                                   outerBounds, "Outer");
  auto outerNode =
      std::make_shared<Experimental::StaticBlueprintNode>(std::move(outerVol));
  root.addChild(outerNode);

  auto trackingGeometry = root.construct({}, gctx, *logger.clone("Geometry"));
  if (trackingGeometry == nullptr) {
    throw std::runtime_error(
        "Failed to build Gen3 tracking geometry in downstream module");
  }

  return trackingGeometry;
}

}  // namespace

ACTS_DEFINE_GEOMETRY_MODULE(buildGeometryModule)
