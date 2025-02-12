// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifierBlueprintNode.hpp"

#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

Volume& GeometryIdentifierBlueprintNode::build(const BlueprintOptions& options,
                                               const GeometryContext& gctx,
                                               const Logger& logger) {
  if (children().size() != 1) {
    throw std::invalid_argument(
        "GeometryIdentifiedBlueprintNode must have exactly one child");
  }

  if (!m_assignment) {
    throw std::invalid_argument(
        "GeometryIdentifiedBlueprintNode has no assignment");
  }

  return children().at(0).build(options, gctx, logger);
}

PortalShellBase& GeometryIdentifierBlueprintNode::connect(
    const BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger) {
  return children().at(0).connect(options, gctx, logger);
}

void GeometryIdentifierBlueprintNode::finalize(const BlueprintOptions& options,
                                               const GeometryContext& gctx,
                                               TrackingVolume& parent,
                                               const Logger& logger) {
  // Run child finalize first!
  children().at(0).finalize(options, gctx, parent, logger);
  for (auto& volume : parent.volumes()) {
    // apply(volume);
    m_assignment(volume);
  }
}

const std::string& GeometryIdentifierBlueprintNode::name() const {
  return m_name;
}

GeometryIdentifierBlueprintNode& GeometryIdentifierBlueprintNode::setLayerTo(
    GeometryIdentifier::Value layer) {
  if (m_assignment) {
    throw std::invalid_argument(
        "GeometryIdentifierBlueprintNode already has an assignment");
  }

  m_assignment = [layer](TrackingVolume& volume) {
    for (auto& surface : volume.surfaces()) {
      GeometryIdentifier id = surface.geometryId();
      id.setLayer(layer);
      surface.assignGeometryId(id);
    }
  };
  return *this;
}

GeometryIdentifierBlueprintNode&
GeometryIdentifierBlueprintNode::incrementLayers() {
  if (m_assignment) {
    throw std::invalid_argument(
        "GeometryIdentifierBlueprintNode already has an assignment");
  }

  struct Incrementer {
    GeometryIdentifier::Value m_value = 0;

    void operator()(TrackingVolume& volume) {
      m_value++;
      GeometryIdentifier id = volume.geometryId();
      id.setLayer(m_value);
      volume.assignGeometryId(id);
    }
  };

  m_assignment = Incrementer();

  return *this;
}

}  // namespace Acts
