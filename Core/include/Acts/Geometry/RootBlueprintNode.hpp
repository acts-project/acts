// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts {

class RootBlueprintNode : public BlueprintNode {
 public:
  struct Config {
    ExtentEnvelope envelope = ExtentEnvelope::Zero();
    GeometryIdentifierHook geometryIdentifierHook = {};
  };

  RootBlueprintNode(const Config& cfg);

  const std::string& name() const override;

  std::unique_ptr<TrackingGeometry> construct(
      const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger());

 protected:
  Volume& build(const Logger& logger = Acts::getDummyLogger()) override;

  PortalShellBase& connect(
      const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  void finalize(TrackingVolume& parent,
                const Logger& logger = Acts::getDummyLogger()) override;

  void addToGraphviz(std::ostream& os) const override;

 private:
  Config m_cfg;
};

}  // namespace Acts
