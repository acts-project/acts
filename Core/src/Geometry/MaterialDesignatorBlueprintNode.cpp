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
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

namespace Acts::Experimental {

namespace {

class CylinderProtoDesignator;
class CuboidProtoDesignator;
class NullDesignator;

class DesignatorBase {
 public:
  virtual ~DesignatorBase() = default;

  virtual void apply(PortalShellBase& shell, const Logger& logger,
                     const std::string& prefix) = 0;

  using FaceVariant =
      std::variant<CylinderVolumeBounds::Face, CuboidVolumeBounds::Face>;

  virtual std::unique_ptr<DesignatorBase> merged(
      const DesignatorBase& other) const = 0;

  virtual std::unique_ptr<DesignatorBase> merged(
      const CylinderProtoDesignator& other) const = 0;

  virtual std::unique_ptr<DesignatorBase> merged(
      const CuboidProtoDesignator& other) const = 0;

  virtual std::unique_ptr<DesignatorBase> merged(
      const NullDesignator& other) const = 0;

  virtual void graphvizLabel(std::ostream& os) const = 0;
};

class NullDesignator : public DesignatorBase {
 public:
  void apply(PortalShellBase& /*shell*/, const Logger& /*logger*/,
             const std::string& /*prefix*/) override {
    throw std::runtime_error("NullDesignator has no apply");
  }

  std::unique_ptr<DesignatorBase> merged(
      const DesignatorBase& other) const override {
    return other.merged(*this);
  }

  std::unique_ptr<DesignatorBase> merged(
      const CylinderProtoDesignator& other) const override;

  std::unique_ptr<DesignatorBase> merged(
      const CuboidProtoDesignator& other) const override;

  std::unique_ptr<DesignatorBase> merged(
      const NullDesignator& /*other*/) const override {
    return std::make_unique<NullDesignator>();
  }

  void graphvizLabel(std::ostream& /*os*/) const override {
    throw std::runtime_error("NullDesignator has no label");
  }
};

class CylinderProtoDesignator : public DesignatorBase {
 public:
  using Face = CylinderVolumeBounds::Face;

  CylinderProtoDesignator(Face face, const DirectedProtoAxis& loc0,
                          const DirectedProtoAxis& loc1,
                          const std::string& prefix) {
    validate(face, loc0, loc1, prefix);

    m_binning.emplace_back(face, loc0, loc1);
  }

  void apply(PortalShellBase& shell, const Logger& logger,
             const std::string& prefix) override {
    auto* cylShell = dynamic_cast<CylinderPortalShell*>(&shell);
    if (cylShell == nullptr) {
      throw std::invalid_argument(prefix +
                                  "Cylinder faces must use a valid face");
    }

    ACTS_DEBUG(prefix << "Binning is set to compatible type");
    using enum CylinderVolumeBounds::Face;

    for (const auto& [face, loc0, loc1] : m_binning) {
      auto* portal = cylShell->portal(face);
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

  void graphvizLabel(std::ostream& os) const override {
    os << "<br/><i>Cylinder Binning</i>";
    for (const auto& [face, loc0, loc1] : m_binning) {
      os << "<br/>" << face;
      os << ": " << loc0.getAxisDirection() << "=" << loc0.getAxis().getNBins();
      os << ", " << loc1.getAxisDirection() << "=" << loc1.getAxis().getNBins();
    }
  }

  std::unique_ptr<DesignatorBase> merged(
      const DesignatorBase& other) const override {
    return other.merged(*this);
  }

  std::unique_ptr<DesignatorBase> merged(
      const CylinderProtoDesignator& other) const override {
    auto merged = std::make_unique<CylinderProtoDesignator>(*this);
    std::ranges::copy(other.m_binning, std::back_inserter(merged->m_binning));
    return merged;
  }

  std::unique_ptr<DesignatorBase> merged(
      const CuboidProtoDesignator& /*other*/) const override {
    throw std::runtime_error(
        "CylinderProtoDesignator cannot be merged with CuboidProtoDesignator");
  }

  std::unique_ptr<DesignatorBase> merged(
      const NullDesignator& /*other*/) const override {
    return std::make_unique<CylinderProtoDesignator>(*this);
  }

 private:
  void validate(Face face, const DirectedProtoAxis& loc0,
                const DirectedProtoAxis& loc1, const std::string& prefix) {
    using enum CylinderVolumeBounds::Face;
    using enum AxisDirection;

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
          throw std::invalid_argument(
              prefix + "Cylinder faces must use (rphi, z) binning");
        }
        break;
      case NegativePhiPlane:
      case PositivePhiPlane:
        throw std::invalid_argument(prefix +
                                    "Phi plane faces are not supported");
        break;
    }
  }

  std::vector<std::tuple<Face, DirectedProtoAxis, DirectedProtoAxis>> m_binning;
};

class CuboidProtoDesignator : public DesignatorBase {
 public:
  using Face = CuboidVolumeBounds::Face;

  CuboidProtoDesignator(Face face, const DirectedProtoAxis& loc0,
                        const DirectedProtoAxis& loc1,
                        const std::string& prefix) {
    validate(loc0, loc1, prefix);

    m_binning.emplace_back(face, loc0, loc1);
  }

  void apply(PortalShellBase& shell, const Logger& logger,
             const std::string& prefix) override {
    auto* cuboidShell = dynamic_cast<CuboidPortalShell*>(&shell);
    if (cuboidShell == nullptr) {
      throw std::invalid_argument(prefix +
                                  "Cuboid faces must use a valid face");
    }

    ACTS_DEBUG(prefix << "Binning is set to compatible type");
    using enum CuboidVolumeBounds::Face;

    for (const auto& [face, loc0, loc1] : m_binning) {
      auto* portal = cuboidShell->portal(face);
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

  void graphvizLabel(std::ostream& os) const override {
    os << "<br/><i>Cuboid Binning</i>";
    for (const auto& [face, loc0, loc1] : m_binning) {
      os << "<br/>" << face;
      os << ": " << loc0.getAxisDirection() << "=" << loc0.getAxis().getNBins();
      os << ", " << loc1.getAxisDirection() << "=" << loc1.getAxis().getNBins();
    }
  }

  std::unique_ptr<DesignatorBase> merged(
      const DesignatorBase& other) const override {
    return other.merged(*this);
  }

  std::unique_ptr<DesignatorBase> merged(
      const CylinderProtoDesignator& /*other*/) const override {
    throw std::runtime_error(
        "CuboidProtoDesignator cannot be merged with CylinderProtoDesignator");
  }

  std::unique_ptr<DesignatorBase> merged(
      const CuboidProtoDesignator& other) const override {
    auto merged = std::make_unique<CuboidProtoDesignator>(*this);
    std::ranges::copy(other.m_binning, std::back_inserter(merged->m_binning));
    return merged;
  }

  std::unique_ptr<DesignatorBase> merged(
      const NullDesignator& /*other*/) const override {
    return std::make_unique<CuboidProtoDesignator>(*this);
  }

 private:
  void validate(const DirectedProtoAxis& loc0, const DirectedProtoAxis& loc1,
                const std::string& prefix) {
    using enum CuboidVolumeBounds::Face;
    using enum AxisDirection;

    // For cuboid faces, the valid axes are always X and Y
    if (loc0.getAxisDirection() != AxisX || loc1.getAxisDirection() != AxisY) {
      throw std::invalid_argument(prefix +
                                  "Cuboid faces must use (x, y) binning");
    }
  }

  std::vector<std::tuple<Face, DirectedProtoAxis, DirectedProtoAxis>> m_binning;
};

std::unique_ptr<DesignatorBase> NullDesignator::merged(
    const CylinderProtoDesignator& other) const {
  return std::make_unique<CylinderProtoDesignator>(other);
}

std::unique_ptr<DesignatorBase> NullDesignator::merged(
    const CuboidProtoDesignator& other) const {
  return std::make_unique<CuboidProtoDesignator>(other);
}

}  // namespace

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
  ss << "" + name() + "";
  ss << "<br/><i>MaterialDesignator</i>";

  impl().m_designator->graphvizLabel(ss);

  os << GraphViz::Node{
      .id = name(), .label = ss.str(), .shape = GraphViz::Shape::Hexagon};
  BlueprintNode::addToGraphviz(os);
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CylinderVolumeBounds::Face face, const DirectedProtoAxis& loc0,
    const DirectedProtoAxis& loc1) {
  impl().m_designator = impl().m_designator->merged(
      CylinderProtoDesignator(face, loc0, loc1, prefix()));

  return *this;
}

MaterialDesignatorBlueprintNode& MaterialDesignatorBlueprintNode::configureFace(
    CuboidVolumeBounds::Face face, const DirectedProtoAxis& loc0,
    const DirectedProtoAxis& loc1) {
  impl().m_designator = impl().m_designator->merged(
      CuboidProtoDesignator(face, loc0, loc1, prefix()));

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
