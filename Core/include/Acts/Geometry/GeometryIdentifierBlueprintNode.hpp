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

#include <memory>

namespace Acts {

class TrackingVolume;
class PortalShellBase;
class Volume;

struct GeometryIdentifierBlueprintNodeImpl;

class GeometryIdentifierBlueprintNode : public BlueprintNode {
 public:
  /// This is not defaulted to it can be pushed to the .cpp file, and we don't
  /// have to define the configuration class.
  ~GeometryIdentifierBlueprintNode() override;
  GeometryIdentifierBlueprintNode();

  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  void finalize(const BlueprintOptions& options, const GeometryContext& gctx,
                TrackingVolume& parent,
                const Logger& logger = Acts::getDummyLogger()) override;

  GeometryIdentifierBlueprintNode& setLayerIdTo(
      GeometryIdentifier::Value layer);
  GeometryIdentifierBlueprintNode& incrementLayerIds(
      GeometryIdentifier::Value start = 0);
  GeometryIdentifierBlueprintNode& setAllVolumeIdsTo(
      GeometryIdentifier::Value volumeId);

  const std::string& name() const override;

 private:
  std::unique_ptr<GeometryIdentifierBlueprintNodeImpl> m_impl;
};

}  // namespace Acts
