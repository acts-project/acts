// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderContainerBlueprintNode.hpp"

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

#include <algorithm>

namespace Acts {

CylinderContainerBlueprintNode::CylinderContainerBlueprintNode(
    const std::string& name, BinningValue direction,
    CylinderVolumeStack::AttachmentStrategy attachmentStrategy,
    CylinderVolumeStack::ResizeStrategy resizeStrategy)
    : m_name(name),
      m_direction(direction),
      m_attachmentStrategy(attachmentStrategy),
      m_resizeStrategy(resizeStrategy) {}

const std::string& CylinderContainerBlueprintNode::name() const {
  return m_name;
}

Volume& CylinderContainerBlueprintNode::build(const Logger& logger) {
  ACTS_DEBUG(prefix() << "cylinder container build");

  if (m_stack.has_value()) {
    ACTS_ERROR(prefix() << "Volume is already built");
    throw std::runtime_error("Volume is already built");
  }

  for (auto& child : children()) {
    Volume& volume = child.build(logger);
    m_childVolumes.push_back(&volume);
    // We need to remember which volume we got from which child, so we can
    // assemble a correct portal shell later
    m_volumeToNode[&volume] = &child;
  }
  ACTS_VERBOSE(prefix() << "-> Collected " << m_childVolumes.size()
                        << " child volumes");

  ACTS_VERBOSE(prefix() << "-> Building the stack");
  m_stack.emplace(m_childVolumes, m_direction, m_attachmentStrategy,
                  m_resizeStrategy, logger);
  ACTS_DEBUG(prefix() << "-> Stack bounds are: " << m_stack->volumeBounds());

  ACTS_DEBUG(prefix() << " *** build complete ***");

  return m_stack.value();
}

CylinderStackPortalShell& CylinderContainerBlueprintNode::connect(
    const GeometryContext& gctx, TrackingVolume& parent, const Logger& logger) {
  ACTS_DEBUG(prefix() << "Cylinder container connect");
  if (!m_stack.has_value()) {
    ACTS_ERROR(prefix() << "Volume is not built");
    throw std::runtime_error("Volume is not built");
  }

  std::vector<CylinderPortalShell*> shells;
  ACTS_VERBOSE(prefix() << "Collecting child shells from " << children().size()
                        << " children");

  // We have child volumes and gaps as bare Volumes in `m_childVolumes` after
  // `build()` has completed. For the stack shell, we need TrackingVolumes in
  // the right order.

  ACTS_VERBOSE(prefix() << "Have " << m_childVolumes.size()
                        << " child volumes");
  for (Volume* volume : m_childVolumes) {
    if (isGapVolume(*volume)) {
      // We need to create a TrackingVolume from the gap and put it in the shell
      auto gapPtr = std::make_unique<TrackingVolume>(*volume);
      ACTS_VERBOSE(prefix() << " - have gap volume: " << gapPtr);
      TrackingVolume& gap = *gapPtr;
      auto& p = m_gapVolumes.emplace_back(
          std::move(gapPtr), std::make_unique<SingleCylinderPortalShell>(gap));

      shells.push_back(p.second.get());
    } else {
      ACTS_VERBOSE(prefix() << "Associate child volume with child node");
      // Figure out which child we got this volume from
      auto it = m_volumeToNode.find(volume);
      if (it == m_volumeToNode.end()) {
        throw std::runtime_error("Volume not found in child volumes");
      }
      BlueprintNode& child = *it->second;

      ACTS_VERBOSE(prefix() << " ~> found child node " << child.name());

      CylinderPortalShell* shell = dynamic_cast<CylinderPortalShell*>(
          &child.connect(gctx, parent, logger));
      if (shell == nullptr) {
        ACTS_ERROR(prefix()
                   << "Child volume of cylinder stack is not a cylinder");
        throw std::runtime_error(
            "Child volume of cylinder stack is not a cylinder");
      }

      shells.push_back(shell);
    }
  }

  // Sanity checks
  throw_assert(shells.size() == m_childVolumes.size(),
               "Number of shells does not match number of child volumes");

  throw_assert(std::ranges::none_of(
                   shells, [](const auto* shell) { return shell == nullptr; }),
               "Invalid shell pointer");

  ACTS_VERBOSE(prefix() << "Producing merged cylinder stack shell in "
                        << m_direction << " direction");
  m_shell.emplace(gctx, std::move(shells), m_direction, logger);

  return m_shell.value();
}

bool CylinderContainerBlueprintNode::isGapVolume(const Volume& volume) const {
  return std::ranges::any_of(
      m_stack->gaps(), [&](const auto& gap) { return gap.get() == &volume; });
}

void CylinderContainerBlueprintNode::visualize(
    IVisualization3D& vis, const GeometryContext& gctx) const {
  if (!m_stack.has_value()) {
    throw std::runtime_error("Cylinder Stack Volume is not built");
  }

  ViewConfig viewConfig{.color = {255, 0, 0}};

  for (const auto& gap : m_stack->gaps()) {
    GeometryView3D::drawVolume(vis, *gap, gctx, Transform3::Identity(),
                               viewConfig);
  }

  BlueprintNode::visualize(vis, gctx);
}

CylinderContainerBlueprintNode& CylinderContainerBlueprintNode::setDirection(
    BinningValue direction) {
  if (m_stack.has_value()) {
    throw std::runtime_error("Cannot change direction after build");
  }
  m_direction = direction;
  return *this;
}

CylinderContainerBlueprintNode&
CylinderContainerBlueprintNode::setAttachmentStrategy(
    CylinderVolumeStack::AttachmentStrategy attachmentStrategy) {
  if (m_stack.has_value()) {
    throw std::runtime_error("Cannot change direction after build");
  }
  m_attachmentStrategy = attachmentStrategy;
  return *this;
}

CylinderContainerBlueprintNode&
CylinderContainerBlueprintNode::setResizeStrategy(
    CylinderVolumeStack::ResizeStrategy resizeStrategy) {
  if (m_stack.has_value()) {
    throw std::runtime_error("Cannot change direction after build");
  }
  m_resizeStrategy = resizeStrategy;
  return *this;
}

void CylinderContainerBlueprintNode::addToGraphviz(std::ostream& os) const {
  GraphViz::Node node{.id = name(),
                      .label = "<b>" + name() + "</b><br/>Cylinder",
                      .shape = GraphViz::Shape::DoubleOctagon};
  os << node << std::endl;
  for (const auto& child : children()) {
    os << indent() << GraphViz::Edge{{.id = name()}, {.id = child.name()}}
       << std::endl;
    child.addToGraphviz(os);
  }
}

}  // namespace Acts
