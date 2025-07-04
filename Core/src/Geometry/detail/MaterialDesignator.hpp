// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"

#include <boost/core/demangle.hpp>

namespace Acts::Experimental::detail {

class CylinderProtoDesignator;
class CuboidProtoDesignator;
class NullDesignator;
template <typename face_enum_t, typename shell_type_t>
class HomogeneousMaterialDesignator;

using CylinderHomogeneousMaterialDesignator =
    HomogeneousMaterialDesignator<CylinderVolumeBounds::Face,
                                  CylinderPortalShell>;

using CuboidHomogeneousMaterialDesignator =
    HomogeneousMaterialDesignator<CuboidVolumeBounds::Face, CuboidPortalShell>;

class DesignatorBase {
 public:
  virtual ~DesignatorBase() = default;

  virtual void apply(PortalShellBase& shell, const Logger& logger,
                     const std::string& prefix) = 0;

  virtual std::string label() const = 0;

  using FaceVariant =
      std::variant<CylinderVolumeBounds::Face, CuboidVolumeBounds::Face>;

  virtual std::unique_ptr<DesignatorBase> merged(
      const DesignatorBase& other) const = 0;

  virtual std::unique_ptr<DesignatorBase> merged(
      const CylinderProtoDesignator& other) const;

  virtual std::unique_ptr<DesignatorBase> merged(
      const CuboidProtoDesignator& other) const;

  virtual std::unique_ptr<DesignatorBase> merged(
      const NullDesignator& other) const;

  virtual std::unique_ptr<DesignatorBase> merged(
      const CylinderHomogeneousMaterialDesignator& other) const;

  virtual std::unique_ptr<DesignatorBase> merged(
      const CuboidHomogeneousMaterialDesignator& other) const;

  virtual void graphvizLabel(std::ostream& os) const = 0;
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

  std::string label() const override { return "CylinderProtoDesignator"; }

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
      os << "<br/> at: " << face;
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

    if (std::ranges::find_if(m_binning, [&](const auto& bin) {
          return std::get<0>(bin) == face;
        }) != m_binning.end()) {
      throw std::invalid_argument(prefix + "Cylinder face already configured");
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
    validate(face, loc0, loc1, prefix);

    m_binning.emplace_back(face, loc0, loc1);
  }

  std::string label() const override { return "CuboidProtoDesignator"; }

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
      os << "<br/> at: " << face;
      os << ": " << loc0.getAxisDirection() << "=" << loc0.getAxis().getNBins();
      os << ", " << loc1.getAxisDirection() << "=" << loc1.getAxis().getNBins();
    }
  }

  std::unique_ptr<DesignatorBase> merged(
      const DesignatorBase& other) const override {
    return other.merged(*this);
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
  void validate(Face face, const DirectedProtoAxis& loc0,
                const DirectedProtoAxis& loc1, const std::string& prefix) {
    using enum CuboidVolumeBounds::Face;
    using enum AxisDirection;

    // For cuboid faces, the valid axes are always X and Y
    if (loc0.getAxisDirection() != AxisX || loc1.getAxisDirection() != AxisY) {
      throw std::invalid_argument(prefix +
                                  "Cuboid faces must use (x, y) binning");
    }

    if (std::ranges::find_if(m_binning, [&](const auto& bin) {
          return std::get<0>(bin) == face;
        }) != m_binning.end()) {
      throw std::invalid_argument(prefix + "Cuboid face already configured");
    }
  }

  std::vector<std::tuple<Face, DirectedProtoAxis, DirectedProtoAxis>> m_binning;
};

template <typename face_enum_t, typename shell_type_t>
class HomogeneousMaterialDesignator : public DesignatorBase {
 public:
  using FaceType = face_enum_t;
  using ShellType = shell_type_t;

  HomogeneousMaterialDesignator(
      FaceType face,
      std::shared_ptr<const HomogeneousSurfaceMaterial> material) {
    m_materials.emplace_back(face, std::move(material));
  }

  std::string label() const override {
    return std::string{shellTypeName()} + " HomogeneousMaterialDesignator";
  }

  void apply(PortalShellBase& shell, const Logger& logger,
             const std::string& prefix) override {
    auto* concreteShell = dynamic_cast<ShellType*>(&shell);
    if (concreteShell == nullptr) {
      ACTS_ERROR(prefix << "Concrete shell type mismatch: configured for "
                        << boost::core::demangle(typeid(ShellType).name())
                        << " but received "
                        << boost::core::demangle(typeid(shell).name()));

      throw std::invalid_argument(prefix + "Concrete shell type mismatch");
    }

    for (const auto& [face, material] : m_materials) {
      auto* portal = concreteShell->portal(face);
      if (portal == nullptr) {
        ACTS_ERROR(prefix << "Portal is nullptr");
        throw std::runtime_error("Portal is nullptr");
      }

      portal->surface().assignSurfaceMaterial(material);
    }
  }

  constexpr std::string_view shellTypeName() const {
    if constexpr (std::is_same_v<ShellType, CylinderPortalShell>) {
      return "Cylinder";
    } else if constexpr (std::is_same_v<ShellType, CuboidPortalShell>) {
      return "Cuboid";
    } else {
      throw std::runtime_error("Unknown shell type");
    }
  }

  void graphvizLabel(std::ostream& os) const override {
    os << "<br/><i>Homog. " << shellTypeName() << " material</i>";
    for (const auto& [face, material] : m_materials) {
      os << "<br/> at: " << face;
    }
  }

  std::unique_ptr<DesignatorBase> merged(
      const DesignatorBase& other) const override {
    return other.merged(*this);
  }

  std::unique_ptr<DesignatorBase> merged(
      const HomogeneousMaterialDesignator& other) const override {
    auto merged = std::make_unique<HomogeneousMaterialDesignator>(*this);
    std::ranges::copy(other.m_materials,
                      std::back_inserter(merged->m_materials));
    return merged;
  }

  std::unique_ptr<DesignatorBase> merged(
      const NullDesignator& /*other*/) const override {
    return std::make_unique<HomogeneousMaterialDesignator>(*this);
  }

 private:
  std::vector<
      std::pair<FaceType, std::shared_ptr<const HomogeneousSurfaceMaterial>>>
      m_materials;
};

class NullDesignator : public DesignatorBase {
 public:
  void apply(PortalShellBase& /*shell*/, const Logger& logger,
             const std::string& prefix) override {
    ACTS_WARNING(prefix << "MaterialDesignator was not configured with any "
                        << "material designation! Check your configuration.");
  }

  std::string label() const override { return "NullDesignator"; }

  std::unique_ptr<DesignatorBase> merged(
      const DesignatorBase& other) const override {
    return other.merged(*this);
  }

  std::unique_ptr<DesignatorBase> merged(
      const CylinderProtoDesignator& other) const override {
    return std::make_unique<CylinderProtoDesignator>(other);
  }

  std::unique_ptr<DesignatorBase> merged(
      const CuboidProtoDesignator& other) const override {
    return std::make_unique<CuboidProtoDesignator>(other);
  }

  std::unique_ptr<DesignatorBase> merged(
      const NullDesignator& /*other*/) const override {
    return std::make_unique<NullDesignator>();
  }

  std::unique_ptr<DesignatorBase> merged(
      const CylinderHomogeneousMaterialDesignator& other) const override {
    return other.merged(*this);
  }

  std::unique_ptr<DesignatorBase> merged(
      const CuboidHomogeneousMaterialDesignator& other) const override {
    return other.merged(*this);
  }

  void graphvizLabel(std::ostream& os) const override {
    os << "<br/><i>NullDesignator</i>";
  }
};

inline std::string mergingError(const DesignatorBase& lhs,
                                const DesignatorBase& rhs) {
  return "MaterialDesignator: Merging of " + lhs.label() + " with " +
         rhs.label() +
         " is not supported. This means you are trying to "
         "configure the same face with different binning or "
         "materials. Please check your configuration.";
}

inline std::unique_ptr<DesignatorBase> DesignatorBase::merged(
    const CylinderProtoDesignator& other) const {
  throw std::runtime_error(mergingError(*this, other));
}

inline std::unique_ptr<DesignatorBase> DesignatorBase::merged(
    const CuboidProtoDesignator& other) const {
  throw std::runtime_error(mergingError(*this, other));
}

inline std::unique_ptr<DesignatorBase> DesignatorBase::merged(
    const NullDesignator& other) const {
  throw std::runtime_error(mergingError(*this, other));
}

inline std::unique_ptr<DesignatorBase> DesignatorBase::merged(
    const CylinderHomogeneousMaterialDesignator& other) const {
  throw std::runtime_error(mergingError(*this, other));
}

inline std::unique_ptr<DesignatorBase> DesignatorBase::merged(
    const CuboidHomogeneousMaterialDesignator& other) const {
  throw std::runtime_error(mergingError(*this, other));
}

}  // namespace Acts::Experimental::detail
