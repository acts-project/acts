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
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include "./detail/MaterialDesignator.hpp"

namespace Acts::Experimental {

namespace detail {
class MaterialDesignatorBlueprintNodeImpl {
 public:
  using CylinderBinning = std::tuple<CylinderVolumeBounds::Face,
                                     DirectedProtoAxis, DirectedProtoAxis>;
  using CuboidBinning = std::tuple<CuboidVolumeBounds::Face, DirectedProtoAxis,
                                   DirectedProtoAxis>;

  using BinningConfig =
      std::variant<std::vector<CylinderBinning>, std::vector<CuboidBinning>>;

  void validateCylinderFaceConfig(CylinderVolumeBounds::Face face,
                                  const DirectedProtoAxis& loc0,
                                  const DirectedProtoAxis& loc1,
                                  const std::string& prefix);

  void validateCuboidFaceConfig(const DirectedProtoAxis& loc0,
                                const DirectedProtoAxis& loc1,
                                const std::string& prefix);

  std::string m_name{};

  std::unique_ptr<DesignatorBase> m_designator{
      std::make_unique<NullDesignator>()};
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

  if (!impl().m_designator) {
    ACTS_ERROR(prefix() << "Designator is not set");
    throw std::runtime_error("Designator is not set");
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

  if (!impl().m_designator) {
    ACTS_ERROR(prefix() << "Designator is not set");
    throw std::runtime_error("Designator is not set");
  }

  auto& shell = children().at(0).connect(options, gctx, logger);

  ACTS_DEBUG(prefix() << "Received shell from child "
                      << children().at(0).name());

  impl().m_designator->apply(shell, logger, prefix());

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
  if (!impl().m_designator) {
    throw std::runtime_error("Binning is not set");
  }

  std::stringstream ss;
  ss << "<b>" + name() + "</b>";
  ss << "<br/>MaterialDesignator";

  impl().m_designator->graphvizLabel(ss);

  os << GraphViz::Node{
      .id = name(), .label = ss.str(), .shape = GraphViz::Shape::Hexagon};
  BlueprintNode::addToGraphviz(os);
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CylinderVolumeBounds::Face face, const DirectedProtoAxis& loc0,
    const DirectedProtoAxis& loc1) {
  impl().m_designator = impl().m_designator->merged(
      detail::CylinderProtoDesignator(face, loc0, loc1, prefix()));

  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CylinderVolumeBounds::Face face,
    std::shared_ptr<const Acts::HomogeneousSurfaceMaterial> material) {
  if (material == nullptr) {
    throw std::invalid_argument(prefix() + "Material is nullptr");
  }

  impl().m_designator = impl().m_designator->merged(
      detail::HomogeneousMaterialDesignator<CylinderVolumeBounds::Face,
                                            CylinderPortalShell>(
          face, std::move(material)));

  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CuboidVolumeBounds::Face face, const DirectedProtoAxis& loc0,
    const DirectedProtoAxis& loc1) {
  impl().m_designator = impl().m_designator->merged(
      detail::CuboidProtoDesignator(face, loc0, loc1, prefix()));

  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CuboidVolumeBounds::Face face,
    std::shared_ptr<const Acts::HomogeneousSurfaceMaterial> material) {
  if (material == nullptr) {
    throw std::invalid_argument(prefix() + "Material is nullptr");
  }

  impl().m_designator = impl().m_designator->merged(
      detail::HomogeneousMaterialDesignator<CuboidVolumeBounds::Face,
                                            CuboidPortalShell>(
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
