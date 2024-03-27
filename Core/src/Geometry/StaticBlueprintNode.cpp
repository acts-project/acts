// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/StaticBlueprintNode.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

namespace Acts {

StaticBlueprintNode::StaticBlueprintNode(const std::string& name,
                                         std::unique_ptr<TrackingVolume> volume)
    : BlueprintNode(name), m_volume(std::move(volume)) {}

Volume& StaticBlueprintNode::build(const Logger& logger) {
  ACTS_DEBUG(prefix() << "static build");
  if (!m_volume) {
    throw std::runtime_error("Volume is not built");
  }

  ACTS_DEBUG(prefix() << "-> returning volume " << *m_volume);
  return *m_volume;
}

void StaticBlueprintNode::connect(TrackingVolume& parent,
                                  const Logger& logger) {
  ACTS_DEBUG(prefix() << "static connect")
  if (!m_volume) {
    throw std::runtime_error("Volume is not built");
  }

  for (auto& child : children()) {
    child.connect(*m_volume);
  }
  parent.addVolume(std::move(m_volume));
}

void StaticBlueprintNode::connect(const Logger& logger) {
  ACTS_DEBUG(prefix() << "static connect")
  if (!m_volume) {
    throw std::runtime_error("Volume is not built");
  }

  for (auto& child : children()) {
    child.connect(*m_volume);
  }
}

void StaticBlueprintNode::visualize(IVisualization3D& vis,
                                    const GeometryContext& gctx) const {
  if (!m_volume) {
    throw std::runtime_error("Volume is not built");
  }

  ViewConfig viewConfig{{100, 100, 100}};

  GeometryView3D::drawVolume(vis, *m_volume, gctx, Transform3::Identity(),
                             viewConfig);
  BlueprintNode::visualize(vis, gctx);
}

}  // namespace Acts
