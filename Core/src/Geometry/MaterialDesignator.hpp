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
#include "Acts/Geometry/DiamondPortalShell.hpp"
#include "Acts/Geometry/DiamondVolumeBounds.hpp"
#include "Acts/Geometry/TrapezoidPortalShell.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"

#include <format>
#include <regex>
#include <string>
#include <variant>

#include <boost/core/demangle.hpp>

namespace Acts::Experimental::detail {

/// Extract the shape name from a portal shell type name, e.g.
/// "Acts::CuboidPortalShell" -> "Cuboid". Falls back to full demangled name if
/// the "PortalShell" suffix is not found.
inline std::string portalShellShapeName(const std::type_info& type) {
  std::string full = boost::core::demangle(type.name());
  std::smatch match;
  static const std::regex regex{R"((\w+)PortalShell$)"};
  if (std::regex_search(full, match, regex)) {
    return match[1].str();
  }
  return full;
}

template <typename face_enum_t, typename shell_type_t>
class ProtoDesignator {
 public:
  using Face = face_enum_t;
  using ShellType = shell_type_t;

  ProtoDesignator(Face face, const DirectedProtoAxis& loc0,
                  const DirectedProtoAxis& loc1, const std::string& prefix) {
    validateAxes(face, loc0, loc1, prefix);
    validateDuplicate(face, prefix);
    m_binning.emplace_back(face, loc0, loc1);
  }

  std::string label() const {
    return std::format("{}ProtoDesignator",
                       portalShellShapeName(typeid(ShellType)));
  }

  void apply(PortalShellBase& shell, const Logger& logger,
             const std::string& prefix) const {
    auto* concreteShell = dynamic_cast<ShellType*>(&shell);
    if (concreteShell == nullptr) {
      throw std::invalid_argument(prefix +
                                  portalShellShapeName(typeid(ShellType)) +
                                  " faces must use a valid face");
    }

    ACTS_DEBUG(prefix << "Binning is set to compatible type");

    for (const auto& [face, loc0, loc1] : m_binning) {
      auto* portal = concreteShell->portal(face).get();
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

  void graphvizLabel(std::ostream& os) const {
    os << "<br/><i>" << portalShellShapeName(typeid(ShellType))
       << " Binning</i>";
    for (const auto& [face, loc0, loc1] : m_binning) {
      os << "<br/> at: " << face;
      os << ": " << loc0.getAxisDirection() << "=" << loc0.getAxis().getNBins();
      os << ", " << loc1.getAxisDirection() << "=" << loc1.getAxis().getNBins();
    }
  }

  ProtoDesignator merged(const ProtoDesignator& other) const {
    ProtoDesignator result = *this;
    std::ranges::copy(other.m_binning, std::back_inserter(result.m_binning));
    return result;
  }

 private:
  void validateAxes(Face face, const DirectedProtoAxis& loc0,
                    const DirectedProtoAxis& loc1, const std::string& prefix) {
    using enum AxisDirection;

    if constexpr (std::is_same_v<ShellType, CylinderPortalShell>) {
      using enum CylinderVolumeBounds::Face;
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
      }
    } else if constexpr (std::is_same_v<ShellType, CuboidPortalShell>) {
      if (loc0.getAxisDirection() != AxisX ||
          loc1.getAxisDirection() != AxisY) {
        throw std::invalid_argument(prefix +
                                    "Cuboid faces must use (x, y) binning");
      }
    } else {
      static_assert(std::is_same_v<ShellType, void>, "Unknown shell type");
    }
  }

  void validateDuplicate(Face face, const std::string& prefix) {
    if (std::ranges::find_if(m_binning, [&](const auto& bin) {
          return std::get<0>(bin) == face;
        }) != m_binning.end()) {
      throw std::invalid_argument(prefix +
                                  portalShellShapeName(typeid(ShellType)) +
                                  " face already configured");
    }
  }

  std::vector<std::tuple<Face, DirectedProtoAxis, DirectedProtoAxis>> m_binning;
};

using CylinderProtoDesignator =
    ProtoDesignator<CylinderVolumeBounds::Face, CylinderPortalShell>;
using CuboidProtoDesignator =
    ProtoDesignator<CuboidVolumeBounds::Face, CuboidPortalShell>;

template <typename face_enum_t, typename shell_type_t>
class ISurfaceMaterialDesignator {
 public:
  using FaceType = face_enum_t;
  using ShellType = shell_type_t;

  ISurfaceMaterialDesignator(FaceType face,
                             std::shared_ptr<const ISurfaceMaterial> material) {
    m_materials.emplace_back(face, std::move(material));
  }

  std::string label() const {
    return std::format("{} ISurfaceMaterialDesignator",
                       portalShellShapeName(typeid(ShellType)));
  }

  void apply(PortalShellBase& shell, const Logger& logger,
             const std::string& prefix) const {
    auto* concreteShell = dynamic_cast<ShellType*>(&shell);
    if (concreteShell == nullptr) {
      ACTS_ERROR(prefix << "Concrete shell type mismatch: configured for "
                        << portalShellShapeName(typeid(ShellType))
                        << " but received "
                        << portalShellShapeName(typeid(shell)));

      throw std::invalid_argument(prefix + "Concrete shell type mismatch");
    }

    for (const auto& [face, material] : m_materials) {
      auto* portal = concreteShell->portal(face).get();
      if (portal == nullptr) {
        ACTS_ERROR(prefix << "Portal is nullptr");
        throw std::runtime_error("Portal is nullptr");
      }

      portal->surface().assignSurfaceMaterial(material);
    }
  }

  void graphvizLabel(std::ostream& os) const {
    os << "<br/><i>Homog. " << portalShellShapeName(typeid(ShellType))
       << " material</i>";
    for (const auto& [face, material] : m_materials) {
      os << "<br/> at: " << face;
    }
  }

  ISurfaceMaterialDesignator merged(
      const ISurfaceMaterialDesignator& other) const {
    ISurfaceMaterialDesignator result = *this;
    std::ranges::copy(other.m_materials,
                      std::back_inserter(result.m_materials));
    return result;
  }

 private:
  std::vector<std::pair<FaceType, std::shared_ptr<const ISurfaceMaterial>>>
      m_materials;
};

using CylinderHomogeneousMaterialDesignator =
    ISurfaceMaterialDesignator<CylinderVolumeBounds::Face, CylinderPortalShell>;

using CuboidHomogeneousMaterialDesignator =
    ISurfaceMaterialDesignator<CuboidVolumeBounds::Face, CuboidPortalShell>;

using TrapezoidHomogeneousMaterialDesignator =
    ISurfaceMaterialDesignator<TrapezoidVolumeBounds::Face,
                               TrapezoidPortalShell>;

using DiamondHomogeneousMaterialDesignator =
    ISurfaceMaterialDesignator<DiamondVolumeBounds::Face, DiamondPortalShell>;

/// Type-erased material designator. std::monostate represents the unconfigured
/// (null) state.
using Designator =
    std::variant<std::monostate, CylinderProtoDesignator, CuboidProtoDesignator,
                 CylinderHomogeneousMaterialDesignator,
                 CuboidHomogeneousMaterialDesignator,
                 TrapezoidHomogeneousMaterialDesignator,
                 DiamondHomogeneousMaterialDesignator>;

/// Merge two designators. std::monostate acts as the identity element.
/// Throws if the two active designator types are incompatible (e.g. cylinder
/// proto vs cuboid proto, or proto vs homogeneous).
inline Designator merge(const Designator& a, const Designator& b) {
  return std::visit(
      [](const auto& x, const auto& y) -> Designator {
        using X = std::decay_t<decltype(x)>;
        using Y = std::decay_t<decltype(y)>;

        if constexpr (std::is_same_v<X, std::monostate>) {
          return y;
        } else if constexpr (std::is_same_v<Y, std::monostate>) {
          return x;
        } else if constexpr (std::is_same_v<X, Y>) {
          return x.merged(y);
        } else {
          throw std::runtime_error(
              std::format("MaterialDesignator: Merging of {} with {} is not "
                          "supported. This means you are trying to configure "
                          "the same face with different binning or materials. "
                          "Please check your configuration.",
                          x.label(), y.label()));
        }
      },
      a, b);
}

/// Apply a designator to a portal shell.
inline void apply(const Designator& d, PortalShellBase& shell,
                  const Logger& logger, const std::string& prefix) {
  std::visit(
      [&](const auto& des) {
        using T = std::decay_t<decltype(des)>;
        if constexpr (std::is_same_v<T, std::monostate>) {
          ACTS_WARNING(prefix << "MaterialDesignator was not configured with "
                              << "any material designation! Check your "
                              << "configuration.");
        } else {
          des.apply(shell, logger, prefix);
        }
      },
      d);
}

/// Write a Graphviz label for a designator to an output stream.
inline void graphvizLabel(const Designator& d, std::ostream& os) {
  std::visit(
      [&](const auto& des) {
        using T = std::decay_t<decltype(des)>;
        if constexpr (std::is_same_v<T, std::monostate>) {
          os << "<br/><i>NullDesignator</i>";
        } else {
          des.graphvizLabel(os);
        }
      },
      d);
}

}  // namespace Acts::Experimental::detail
