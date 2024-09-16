// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/StaticBlueprintNode.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

namespace Acts {

StaticBlueprintNode::StaticBlueprintNode(std::unique_ptr<TrackingVolume> volume)
    : m_volume(std::move(volume)) {}

Volume& StaticBlueprintNode::build(const Logger& logger) {
  ACTS_DEBUG(prefix() << "static build");
  if (!m_volume) {
    throw std::runtime_error("Volume is not built");
  }

  ACTS_DEBUG(prefix() << "Building volume (" << name() << ") with "
                      << children().size() << " children");
  for (auto& child : children()) {
    child.build(logger);
  }

  ACTS_DEBUG(prefix() << "-> returning volume " << *m_volume);
  return *m_volume;
}

PortalShellBase& StaticBlueprintNode::connect(const GeometryContext& gctx,
                                              TrackingVolume& parent,
                                              const Logger& logger) {
  ACTS_DEBUG(prefix() << "Static connect");
  if (m_volume == nullptr) {
    throw std::runtime_error("Volume is not present");
  }

  ACTS_DEBUG(prefix() << "Connecting parent volume (" << name() << ") with "
                      << children().size() << " children");

  for (auto& child : children()) {
    auto& shell = child.connect(gctx, parent, logger);
    // Register ourselves on the outside of the shell
    shell.connectOuter(*m_volume);
  }

  VolumeBounds::BoundsType type = m_volume->volumeBounds().type();
  if (type == VolumeBounds::eCylinder) {
    m_shell = std::make_unique<SingleCylinderPortalShell>(*m_volume);

  } else if (type == VolumeBounds::eCuboid) {
    throw std::logic_error("Cuboid is not implemented yet");

  } else {
    throw std::logic_error("Volume type is not supported");
  }

  ACTS_DEBUG(prefix() << " Adding volume (" << m_volume->volumeName()
                      << ") to parent volume (" << parent.volumeName() << ")");
  parent.addVolume(std::move(m_volume));

  assert(m_shell != nullptr);
  return *m_shell;
}

void StaticBlueprintNode::visualize(IVisualization3D& vis,
                                    const GeometryContext& gctx) const {
  if (!m_volume) {
    throw std::runtime_error("Volume is not built");
  }

  ViewConfig viewConfig{.color = {100, 100, 100}};

  GeometryView3D::drawVolume(vis, *m_volume, gctx, Transform3::Identity(),
                             viewConfig);
  BlueprintNode::visualize(vis, gctx);
}

const std::string& StaticBlueprintNode::name() const {
  static const std::string uninitialized = "uninitialized";
  if (m_volume == nullptr) {
    return uninitialized;
  }
  return m_volume->volumeName();
}

void StaticBlueprintNode::addToGraphviz(std::ostream& os) const {
  std::stringstream ss;
  ss << "<b>" << name() << "</b>";
  ss << "<br/>";
  switch (m_volume->volumeBounds().type()) {
    case VolumeBounds::eCylinder:
      ss << "Cylinder";
      break;
    case VolumeBounds::eCuboid:
      ss << "Cuboid";
      break;
    case VolumeBounds::eCone:
      ss << "Cone";
      break;
    case VolumeBounds::eCutoutCylinder:
      ss << "CutoutCylinder";
      break;
    case VolumeBounds::eGenericCuboid:
      ss << "GenericCuboid";
      break;
    case VolumeBounds::eTrapezoid:
      ss << "Trapezoid";
      break;
    default:
      ss << "Other";
  }

  GraphViz::Node node{
      .id = name(), .label = ss.str(), .shape = GraphViz::Shape::Rectangle};

  os << node;

  BlueprintNode::addToGraphviz(os);
}

}  // namespace Acts
