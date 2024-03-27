// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderContainerBlueprintNode.hpp"

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

namespace Acts {

CylinderContainerBlueprintNode::CylinderContainerBlueprintNode(
    const std::string& name, BinningValue direction,
    CylinderVolumeStack::AttachmentStrategy attachmentStrategy,
    CylinderVolumeStack::ResizeStrategy resizeStrategy)
    : BlueprintNode(name),
      m_direction(direction),
      m_attachmentStrategy(attachmentStrategy),
      m_resizeStrategy(resizeStrategy) {}

Volume& CylinderContainerBlueprintNode::build(const Logger& logger) {
  ACTS_DEBUG(prefix() << "cylinder container build");

  if (m_stack.has_value()) {
    throw std::runtime_error("Volume is already built");
  }

  for (auto& child : children()) {
    m_childVolumes.push_back(&child.build(logger));
  }
  ACTS_VERBOSE(prefix() << "-> Collected " << m_childVolumes.size()
                        << " child volumes");

  ACTS_VERBOSE(prefix() << "-> Building the stack");
  m_stack.emplace(m_childVolumes, m_direction, m_attachmentStrategy,
                  m_resizeStrategy, logger);
  ACTS_DEBUG(prefix() << "-> Stack bounds are: " << m_stack->volumeBounds());

  return m_stack.value();
}

void CylinderContainerBlueprintNode::connect(TrackingVolume& parent,
                                             const Logger& logger) {
  ACTS_DEBUG(prefix() << "cylinder container connect");
  if (!m_stack.has_value()) {
    throw std::runtime_error("Volume is not built");
  }

  for (auto& child : children()) {
    child.connect(parent);
  }

  for (auto& gap : m_stack->gaps()) {
    auto tv = std::make_unique<TrackingVolume>(*gap);
    parent.addVolume(std::move(tv));
  }
}

void CylinderContainerBlueprintNode::visualize(
    IVisualization3D& vis, const GeometryContext& gctx) const {
  if (!m_stack.has_value()) {
    throw std::runtime_error("Cylinder Stack Volume is not built");
  }

  ViewConfig viewConfig{{255, 0, 0}};

  for (const auto& gap : m_stack->gaps()) {
    GeometryView3D::drawVolume(vis, *gap, gctx, Transform3::Identity(),
                               viewConfig);
  }

  BlueprintNode::visualize(vis, gctx);
}

}  // namespace Acts
