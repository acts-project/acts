// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"

#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DiamondPortalShell.hpp"
#include "Acts/Geometry/DiamondVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/TrapezoidPortalShell.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include "./MaterialDesignator.hpp"

namespace Acts::Experimental {

namespace detail {
class MaterialDesignatorBlueprintNodeImpl {
 public:
  std::string m_name{};
  Designator m_designator{std::monostate{}};
};

}  // namespace detail

MaterialDesignatorBlueprintNode::MaterialDesignatorBlueprintNode(
    const std::string& name) {
  m_impl = std::make_unique<detail::MaterialDesignatorBlueprintNodeImpl>();
  m_impl->m_name = name;
}

const std::string& MaterialDesignatorBlueprintNode::name() const {
  return impl().m_name;
}

void MaterialDesignatorBlueprintNode::toStream(std::ostream& os) const {
  os << "MaterialDesignatorBlueprintNode(" << name() << ")";
}

Volume& MaterialDesignatorBlueprintNode::build(const BlueprintOptions& options,
                                               const GeometryContext& gctx,
                                               const Logger& logger) {
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "MaterialDesignatorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "MaterialDesignatorBlueprintNode must have exactly one child");
  }

  return children().at(0).build(options, gctx, logger);
}

PortalShellBase& MaterialDesignatorBlueprintNode::connect(
    const BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger) {
  ACTS_DEBUG(prefix() << "MaterialDesignatorBlueprintNode::connect");
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "MaterialDesignatorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "MaterialDesignatorBlueprintNode must have exactly one child");
  }

  auto& shell = children().at(0).connect(options, gctx, logger);

  ACTS_DEBUG(prefix() << "Received shell from child "
                      << children().at(0).name());

  detail::apply(impl().m_designator, shell, logger, prefix());

  return shell;
}

void MaterialDesignatorBlueprintNode::finalize(const BlueprintOptions& options,
                                               const GeometryContext& gctx,
                                               TrackingVolume& parent,
                                               const Logger& logger) {
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "MaterialDesignatorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "MaterialDesignatorBlueprintNode must have exactly one child");
  }
  return children().at(0).finalize(options, gctx, parent, logger);
}

void MaterialDesignatorBlueprintNode::addToGraphviz(std::ostream& os) const {
  std::stringstream ss;
  ss << "<b>" + name() + "</b>";
  ss << "<br/>MaterialDesignator";

  detail::graphvizLabel(impl().m_designator, ss);

  os << GraphViz::Node{
      .id = name(), .label = ss.str(), .shape = GraphViz::Shape::Hexagon};
  BlueprintNode::addToGraphviz(os);
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CylinderVolumeBounds::Face face, const DirectedProtoAxis& loc0,
    const DirectedProtoAxis& loc1) {
  impl().m_designator = detail::merge(
      impl().m_designator,
      detail::CylinderProtoDesignator(face, loc0, loc1, prefix()));
  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CylinderVolumeBounds::Face face,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  if (material == nullptr) {
    throw std::invalid_argument(prefix() + "Material is nullptr");
  }
  impl().m_designator = detail::merge(
      impl().m_designator,
      detail::CylinderHomogeneousMaterialDesignator(face, std::move(material)));
  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CuboidVolumeBounds::Face face, const DirectedProtoAxis& loc0,
    const DirectedProtoAxis& loc1) {
  impl().m_designator = detail::merge(
      impl().m_designator,
      detail::CuboidProtoDesignator(face, loc0, loc1, prefix()));
  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CuboidVolumeBounds::Face face,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  if (material == nullptr) {
    throw std::invalid_argument(prefix() + "Material is nullptr");
  }
  impl().m_designator = detail::merge(
      impl().m_designator,
      detail::CuboidHomogeneousMaterialDesignator(face, std::move(material)));
  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    TrapezoidVolumeBounds::Face face,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  if (material == nullptr) {
    throw std::invalid_argument(prefix() + "Material is nullptr");
  }
  impl().m_designator = impl().m_designator->merged(
      detail::ISurfaceMaterialDesignator<TrapezoidVolumeBounds::Face,
                                         TrapezoidPortalShell>(
          face, std::move(material)));
  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    DiamondVolumeBounds::Face face,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  if (material == nullptr) {
    throw std::invalid_argument(prefix() + "Material is nullptr");
  }
  impl().m_designator = impl().m_designator->merged(
      detail::ISurfaceMaterialDesignator<DiamondVolumeBounds::Face,
                                         DiamondPortalShell>(
          face, std::move(material)));
  return *this;
}

detail::MaterialDesignatorBlueprintNodeImpl&
MaterialDesignatorBlueprintNode::impl() {
  if (!m_impl) {
    throw std::runtime_error("MaterialDesignatorBlueprintNodeImpl is not set");
  }
  return *m_impl;
}

const detail::MaterialDesignatorBlueprintNodeImpl&
MaterialDesignatorBlueprintNode::impl() const {
  if (!m_impl) {
    throw std::runtime_error("MaterialDesignatorBlueprintNodeImpl is not set");
  }
  return *m_impl;
}

MaterialDesignatorBlueprintNode::~MaterialDesignatorBlueprintNode() = default;

}  // namespace Acts::Experimental
