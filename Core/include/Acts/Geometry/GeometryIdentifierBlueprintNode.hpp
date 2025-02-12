// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class TrackingVolume;
class PortalShellBase;
class Volume;

class GeometryIdentifierBlueprintNode : public BlueprintNode {
 public:
  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  void finalize(const BlueprintOptions& options, const GeometryContext& gctx,
                TrackingVolume& parent,
                const Logger& logger = Acts::getDummyLogger()) override;

  using AssignmentFunction = std::function<void(TrackingVolume& volume)>;

  //   GeometryIdentifierBlueprintNode& volume(GeometryIdentifier::Volume
  //   volume);
  GeometryIdentifierBlueprintNode& setLayerTo(GeometryIdentifier::Value layer);
  GeometryIdentifierBlueprintNode& incrementLayers();

  const std::string& name() const override;

 protected:
  AssignmentFunction m_assignment;

  std::string m_name = "DummyGeoId";
};

}  // namespace Acts
