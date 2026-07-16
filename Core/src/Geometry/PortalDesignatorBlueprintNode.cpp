// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/PortalDesignatorBlueprintNode.hpp"

#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Utilities/GraphViz.hpp"

#include "./PortalDesignator.hpp"

namespace Acts::Experimental {

namespace detail {
class PortalDesignatorBlueprintNodeImpl {
 public:
  std::string m_name{};
  // Monostate is the null state, where no tags are configured.
  PortalTagDesignatorVariant m_designator{std::monostate{}};
  // The child shell captured during connect(), tagged during finalize().
  PortalShellBase* m_shell{nullptr};
};

}  // namespace detail

PortalDesignatorBlueprintNode::PortalDesignatorBlueprintNode(
    std::string_view name) {
  m_impl = std::make_unique<detail::PortalDesignatorBlueprintNodeImpl>();
  m_impl->m_name = name;
}

const std::string& PortalDesignatorBlueprintNode::name() const {
  return impl().m_name;
}

void PortalDesignatorBlueprintNode::toStream(std::ostream& os) const {
  os << "PortalDesignatorBlueprintNode(" << name() << ")";
}

Volume& PortalDesignatorBlueprintNode::build(const BlueprintOptions& options,
                                             const GeometryContext& gctx,
                                             const Logger& logger) {
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "PortalDesignatorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "PortalDesignatorBlueprintNode must have exactly one child");
  }

  return children().at(0).build(options, gctx, logger);
}

PortalShellBase& PortalDesignatorBlueprintNode::connect(
    const BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger) {
  ACTS_DEBUG(prefix() << "PortalDesignatorBlueprintNode::connect");
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "PortalDesignatorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "PortalDesignatorBlueprintNode must have exactly one child");
  }

  auto& shell = children().at(0).connect(options, gctx, logger);

  ACTS_DEBUG(prefix() << "Received shell from child "
                      << children().at(0).name());

  // Capture the shell so we can apply tags in finalize(). Tagging cannot happen
  // here: parent container nodes may still merge/fuse this shell's portals,
  // which replaces the portal objects. By finalize(), the stacking has written
  // the final shared portals back into the shell's face slots.
  impl().m_shell = &shell;

  return shell;
}

void PortalDesignatorBlueprintNode::finalize(const BlueprintOptions& options,
                                             const GeometryContext& gctx,
                                             TrackingVolume& parent,
                                             const Logger& logger) {
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "PortalDesignatorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "PortalDesignatorBlueprintNode must have exactly one child");
  }

  if (impl().m_shell == nullptr) {
    ACTS_ERROR(prefix() << "No shell captured during connect()");
    throw std::runtime_error(
        "PortalDesignatorBlueprintNode: no shell captured during connect()");
  }

  detail::applyTags(impl().m_designator, *impl().m_shell, logger, prefix());

  children().at(0).finalize(options, gctx, parent, logger);
}

PortalDesignatorBlueprintNode& PortalDesignatorBlueprintNode::tagFace(
    CylinderVolumeBounds::Face face, const std::string& label) {
  impl().m_designator = detail::mergeTags(
      impl().m_designator,
      detail::CylinderPortalTagDesignator(face, label, prefix()));
  return *this;
}

PortalDesignatorBlueprintNode& PortalDesignatorBlueprintNode::tagFace(
    CuboidVolumeBounds::Face face, const std::string& label) {
  impl().m_designator = detail::mergeTags(
      impl().m_designator,
      detail::CuboidPortalTagDesignator(face, label, prefix()));
  return *this;
}

void PortalDesignatorBlueprintNode::addToGraphviz(std::ostream& os) const {
  std::stringstream ss;
  ss << "<b>" + name() + "</b>";
  ss << "<br/>PortalDesignator";

  detail::graphvizLabelTags(impl().m_designator, ss);

  os << GraphViz::Node{
      .id = name(), .label = ss.str(), .shape = GraphViz::Shape::Hexagon};
  BlueprintNode::addToGraphviz(os);
}

detail::PortalDesignatorBlueprintNodeImpl&
PortalDesignatorBlueprintNode::impl() {
  if (!m_impl) {
    throw std::runtime_error("PortalDesignatorBlueprintNodeImpl is not set");
  }
  return *m_impl;
}

const detail::PortalDesignatorBlueprintNodeImpl&
PortalDesignatorBlueprintNode::impl() const {
  if (!m_impl) {
    throw std::runtime_error("PortalDesignatorBlueprintNodeImpl is not set");
  }
  return *m_impl;
}

PortalDesignatorBlueprintNode::~PortalDesignatorBlueprintNode() = default;

}  // namespace Acts::Experimental
