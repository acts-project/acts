// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

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

  explicit MaterialDesignatorBlueprintNodeImpl(
      MaterialDesignatorBlueprintNode& parent)
      : m_parent(&parent) {}

  void handleCylinderBinning(
      CylinderPortalShell& cylShell,
      const std::vector<std::tuple<CylinderPortalShell::Face, DirectedProtoAxis,
                                   DirectedProtoAxis>>& binning,
      const Logger& logger, const std::string& prefix);

  void handleCuboidBinning(
      CuboidPortalShell& cuboidShell,
      const std::vector<std::tuple<CuboidVolumeBounds::Face, DirectedProtoAxis,
                                   DirectedProtoAxis>>& binning,
      const Logger& logger, const std::string& prefix);

  void validateCylinderFaceConfig(CylinderVolumeBounds::Face face,
                                  const DirectedProtoAxis& loc0,
                                  const DirectedProtoAxis& loc1,
                                  const std::string& prefix);

  void validateCuboidFaceConfig(const DirectedProtoAxis& loc0,
                                const DirectedProtoAxis& loc1,
                                const std::string& prefix);

  PortalShellBase& connect(const BlueprintOptions& options,
                           const GeometryContext& gctx, const Logger& logger,
                           const std::string& prefix);

  auto children() const {
    const auto& p = *m_parent;
    return p.children();
  }

  auto children() { return m_parent->children(); }

  MaterialDesignatorBlueprintNode* m_parent{nullptr};

  std::string m_name{};

  std::optional<BinningConfig> m_binning{};

  std::unique_ptr<detail::MaterialDesignatorBlueprintNodeImpl> m_impl;
};

void MaterialDesignatorBlueprintNodeImpl::validateCylinderFaceConfig(
    CylinderVolumeBounds::Face face, const DirectedProtoAxis& loc0,
    const DirectedProtoAxis& loc1, const std::string& prefix) {
  using enum CylinderVolumeBounds::Face;
  using enum AxisDirection;

  // Check if we already have a different volume type configured
  if (m_binning.has_value() &&
      !std::holds_alternative<std::vector<std::tuple<
          CylinderVolumeBounds::Face, DirectedProtoAxis, DirectedProtoAxis>>>(
          m_binning.value())) {
    throw std::invalid_argument(
        prefix + "Cannot mix volume types in material configuration");
  }

  // Validate axis directions based on face type
  switch (face) {
    case NegativeDisc:
    case PositiveDisc:
      if (loc0.getAxisDirection() != AxisR ||
          loc1.getAxisDirection() != AxisPhi) {
        throw std::invalid_argument(prefix +
                                    "Disc faces must use (r, phi) binning");
      }
      break;
    case OuterCylinder:
    case InnerCylinder:
      if (loc0.getAxisDirection() != AxisRPhi ||
          loc1.getAxisDirection() != AxisZ) {
        throw std::invalid_argument(prefix +
                                    "Cylinder faces must use (phi, z) binning");
      }
      break;
    case NegativePhiPlane:
    case PositivePhiPlane:
      throw std::invalid_argument(prefix + "Phi plane faces are not supported");
      break;
  }
}

void MaterialDesignatorBlueprintNodeImpl::handleCylinderBinning(
    CylinderPortalShell& cylShell,
    const std::vector<std::tuple<CylinderVolumeBounds::Face, DirectedProtoAxis,
                                 DirectedProtoAxis>>& binning,
    const Logger& logger, const std::string& prefix) {
  ACTS_DEBUG(prefix << "Binning is set to compatible type");
  using enum CylinderVolumeBounds::Face;

  for (const auto& [face, loc0, loc1] : binning) {
    auto* portal = cylShell.portal(face);
    if (portal == nullptr) {
      ACTS_ERROR(prefix << "Portal is nullptr");
      throw std::runtime_error("Portal is nullptr");
    }

    ACTS_DEBUG(prefix << "Assigning material with binning: " << loc0 << ", "
                      << loc1 << " to face " << face);

    portal->surface().assignSurfaceMaterial(
        std::make_shared<ProtoGridSurfaceMaterial>(std::vector{loc0, loc1}));
  }
}

void MaterialDesignatorBlueprintNodeImpl::validateCuboidFaceConfig(
    const DirectedProtoAxis& loc0, const DirectedProtoAxis& loc1,
    const std::string& prefix) {
  using enum CuboidVolumeBounds::Face;
  using enum AxisDirection;

  // Check if we already have a different volume type configured
  if (m_binning.has_value() &&
      !std::holds_alternative<std::vector<std::tuple<
          CuboidVolumeBounds::Face, DirectedProtoAxis, DirectedProtoAxis>>>(
          m_binning.value())) {
    throw std::invalid_argument(
        prefix + "Cannot mix volume types in material configuration");
  }

  // For cuboid faces, the valid axes are always X and Y
  if (loc0.getAxisDirection() != AxisX || loc1.getAxisDirection() != AxisY) {
    throw std::invalid_argument(prefix +
                                "Cuboid faces must use (x, y) binning");
  }
}

void MaterialDesignatorBlueprintNodeImpl::handleCuboidBinning(
    CuboidPortalShell& cuboidShell,
    const std::vector<std::tuple<CuboidVolumeBounds::Face, DirectedProtoAxis,
                                 DirectedProtoAxis>>& binning,
    const Logger& logger, const std::string& prefix) {
  ACTS_DEBUG(prefix << "Binning is set to compatible type");
  using enum CuboidVolumeBounds::Face;

  for (const auto& [face, loc0, loc1] : binning) {
    auto* portal = cuboidShell.portal(face);
    if (portal == nullptr) {
      ACTS_ERROR(prefix << "Portal is nullptr");
      throw std::runtime_error("Portal is nullptr");
    }

    ACTS_DEBUG(prefix << "Assigning material with binning: " << loc0 << ", "
                      << loc1 << " to face " << face);

    portal->surface().assignSurfaceMaterial(
        std::make_shared<ProtoGridSurfaceMaterial>(std::vector{loc0, loc1}));
  }
}

PortalShellBase& MaterialDesignatorBlueprintNodeImpl::connect(
    const BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger, const std::string& prefix) {
  ACTS_DEBUG(prefix << "MaterialDesignatorBlueprintNode::connect");
  if (children().size() != 1) {
    ACTS_ERROR(prefix << "MaterialDesignatorBlueprintNode must have exactly "
                         "one child, but has "
                      << children().size());
    throw std::runtime_error(
        "MaterialDesignatorBlueprintNode must have exactly one child");
  }
  if (!m_binning) {
    ACTS_ERROR(prefix << "Binning is not set");
    throw std::runtime_error("Binning is not set");
  }

  auto& shell = children().at(0).connect(options, gctx, logger);

  ACTS_DEBUG(prefix << "Received shell from child " << children().at(0).name());

  if (!m_binning.has_value()) {
    ACTS_WARNING(
        prefix << "Material designator node has no assignment configured");
    return shell;
  }

  if (auto* cylShell = dynamic_cast<CylinderPortalShell*>(&shell)) {
    if (auto* cylBinning = std::get_if<std::vector<std::tuple<
            CylinderVolumeBounds::Face, DirectedProtoAxis, DirectedProtoAxis>>>(
            &m_binning.value());
        cylBinning != nullptr) {
      handleCylinderBinning(*cylShell, *cylBinning, logger, prefix);
    } else {
      ACTS_ERROR(prefix << "Binning is set to unknown type");
      throw std::runtime_error("Unknown binning type");
    }
  } else if (auto* cuboidShell = dynamic_cast<CuboidPortalShell*>(&shell)) {
    if (auto* cuboidBinning = std::get_if<std::vector<std::tuple<
            CuboidVolumeBounds::Face, DirectedProtoAxis, DirectedProtoAxis>>>(
            &m_binning.value());
        cuboidBinning != nullptr) {
      handleCuboidBinning(*cuboidShell, *cuboidBinning, logger, prefix);
    } else {
      ACTS_ERROR(prefix << "Binning is set to unknown type");
      throw std::runtime_error("Unknown binning type");
    }
  } else {
    ACTS_ERROR(prefix << "Shell is not supported");
    throw std::runtime_error("Shell is not supported");
  }

  return shell;
}

}  // namespace detail

MaterialDesignatorBlueprintNode::MaterialDesignatorBlueprintNode(
    const std::string& name) {
  m_impl = std::make_unique<detail::MaterialDesignatorBlueprintNodeImpl>(*this);
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

  if (!impl().m_binning) {
    ACTS_ERROR(prefix() << "Binning is not set");
    throw std::runtime_error("Binning is not set");
  }

  return children().at(0).build(options, gctx, logger);
}

PortalShellBase& MaterialDesignatorBlueprintNode::connect(
    const BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger) {
  return impl().connect(options, gctx, logger, prefix());
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
  if (!impl().m_binning) {
    throw std::runtime_error("Binning is not set");
  }

  std::stringstream ss;
  ss << "" + name() + "";
  ss << "<br/><i>MaterialDesignator</i>";

  std::visit(
      overloaded{[&](const std::vector<
                     std::tuple<CylinderPortalShell::Face, DirectedProtoAxis,
                                DirectedProtoAxis>>& binning) {
                   ss << "<br/><i>Cylinder Binning</i>";
                   for (const auto& [face, loc0, loc1] : binning) {
                     ss << "<br/>" << face;
                     ss << ": " << loc0.getAxisDirection() << "="
                        << loc0.getAxis().getNBins();
                     ss << ", " << loc1.getAxisDirection() << "="
                        << loc1.getAxis().getNBins();
                   }
                 },
                 [&](const std::vector<
                     std::tuple<CuboidVolumeBounds::Face, DirectedProtoAxis,
                                DirectedProtoAxis>>& binning) {
                   ss << "<br/><i>Cuboid Binning</i>";
                   for (const auto& [face, loc0, loc1] : binning) {
                     ss << "<br/>" << face;
                     ss << ": " << loc0.getAxisDirection() << "="
                        << loc0.getAxis().getNBins();
                     ss << ", " << loc1.getAxisDirection() << "="
                        << loc1.getAxis().getNBins();
                   }
                 },
                 [](const auto& /*binning*/) {
                   // No output in all other cases
                 }},
      impl().m_binning.value());
  os << GraphViz::Node{
      .id = name(), .label = ss.str(), .shape = GraphViz::Shape::Hexagon};
  BlueprintNode::addToGraphviz(os);
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CylinderVolumeBounds::Face face, const DirectedProtoAxis& loc0,
    const DirectedProtoAxis& loc1) {
  impl().validateCylinderFaceConfig(face, loc0, loc1, prefix());

  if (!impl().m_binning) {
    impl().m_binning =
        std::vector<std::tuple<CylinderVolumeBounds::Face, DirectedProtoAxis,
                               DirectedProtoAxis>>{};
  }

  auto& binning =
      std::get<std::vector<std::tuple<CylinderVolumeBounds::Face,
                                      DirectedProtoAxis, DirectedProtoAxis>>>(
          impl().m_binning.value());
  binning.emplace_back(face, loc0, loc1);

  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CuboidVolumeBounds::Face face, const DirectedProtoAxis& loc0,
    const DirectedProtoAxis& loc1) {
  impl().validateCuboidFaceConfig(loc0, loc1, prefix());

  if (!impl().m_binning) {
    impl().m_binning =
        std::vector<std::tuple<CuboidVolumeBounds::Face, DirectedProtoAxis,
                               DirectedProtoAxis>>{};
  }

  auto& binning =
      std::get<std::vector<std::tuple<CuboidVolumeBounds::Face,
                                      DirectedProtoAxis, DirectedProtoAxis>>>(
          impl().m_binning.value());
  binning.emplace_back(face, loc0, loc1);

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
