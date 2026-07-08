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
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <format>
#include <ostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <variant>
#include <vector>

#include <boost/core/demangle.hpp>

namespace Acts::Experimental::detail {

/// Extract the shape name from a portal shell type name, e.g.
/// "Acts::CuboidPortalShell" -> "Cuboid". Falls back to full demangled name if
/// the "PortalShell" suffix is not found.
inline std::string portalShellShapeName(const std::type_info& type) {
  static const std::regex kPortalShellRegex{R"((\w+)PortalShell$)"};
  std::string full = boost::core::demangle(type.name());
  if (std::smatch match; std::regex_search(full, match, kPortalShellRegex)) {
    return match[1].str();
  }
  return full;
}

/// Designator that assigns string tags to specific faces of a portal shell.
///
/// Unlike the material designator, tagging happens in the *finalize* phase of
/// the blueprint construction (after all portal merging/fusing), so the portal
/// retrieved here is the final, shared portal that ends up in the geometry.
/// @tparam face_enum_t The volume-bounds face enum (e.g. CylinderVolumeBounds::Face)
/// @tparam shell_type_t The concrete portal shell base type (e.g. CylinderPortalShell)
template <typename face_enum_t, typename shell_type_t>
class PortalTagDesignator {
 public:
  using Face = face_enum_t;
  using ShellType = shell_type_t;

  PortalTagDesignator(Face face, std::string label, const std::string& prefix) {
    validateDuplicate(face, prefix);
    m_tags.emplace_back(face, std::move(label));
  }

  std::string label() const {
    return std::format("{}PortalTagDesignator",
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

    for (const auto& [face, tag] : m_tags) {
      auto* portal = concreteShell->portal(face).get();
      if (portal == nullptr) {
        ACTS_ERROR(prefix << "Portal is nullptr");
        throw std::runtime_error("Portal is nullptr");
      }

      ACTS_DEBUG(prefix << "Tagging face " << face << " with '" << tag << "'");
      portal->addTag(tag);
    }
  }

  void graphvizLabel(std::ostream& os) const {
    os << "<br/><i>" << portalShellShapeName(typeid(ShellType)) << " tags</i>";
    for (const auto& [face, tag] : m_tags) {
      os << "<br/> at: " << face << ": " << tag;
    }
  }

  PortalTagDesignator merged(const PortalTagDesignator& other) const {
    PortalTagDesignator result = *this;
    std::ranges::copy(other.m_tags, std::back_inserter(result.m_tags));
    return result;
  }

 private:
  void validateDuplicate(Face face, const std::string& prefix) {
    if (std::ranges::find_if(m_tags, [&](const auto& entry) {
          return entry.first == face;
        }) != m_tags.end()) {
      throw std::invalid_argument(prefix +
                                  portalShellShapeName(typeid(ShellType)) +
                                  " face already tagged");
    }
  }

  std::vector<std::pair<Face, std::string>> m_tags;
};

using CylinderPortalTagDesignator =
    PortalTagDesignator<CylinderVolumeBounds::Face, CylinderPortalShell>;
using CuboidPortalTagDesignator =
    PortalTagDesignator<CuboidVolumeBounds::Face, CuboidPortalShell>;

/// Type-erased portal tag designator. std::monostate represents the
/// unconfigured (null) state.
using PortalTagDesignatorVariant =
    std::variant<std::monostate, CylinderPortalTagDesignator,
                 CuboidPortalTagDesignator>;

/// Merge two designators. std::monostate acts as the identity element.
/// Throws if the two active designator types are incompatible (e.g. cylinder
/// vs cuboid).
inline PortalTagDesignatorVariant mergeTags(
    const PortalTagDesignatorVariant& a, const PortalTagDesignatorVariant& b) {
  return std::visit(
      []<typename X, typename Y>(const X& x,
                                 const Y& y) -> PortalTagDesignatorVariant {
        if constexpr (std::is_same_v<X, std::monostate>) {
          return y;
        } else if constexpr (std::is_same_v<Y, std::monostate>) {
          return x;
        } else if constexpr (std::is_same_v<X, Y>) {
          return x.merged(y);
        } else {
          throw std::runtime_error(std::format(
              "PortalDesignator: Mixing {} with {} is not supported. A portal "
              "designator node can only tag faces of a single volume shape.",
              x.label(), y.label()));
        }
      },
      a, b);
}

/// Apply a designator to a portal shell.
inline void applyTags(const PortalTagDesignatorVariant& d,
                      PortalShellBase& shell, const Logger& logger,
                      const std::string& prefix) {
  std::visit(
      [&]<typename T>(const T& des) {
        if constexpr (std::is_same_v<T, std::monostate>) {
          ACTS_WARNING(prefix << "PortalDesignator was not configured with any "
                              << "face tags! Check your configuration.");
        } else {
          des.apply(shell, logger, prefix);
        }
      },
      d);
}

/// Write a Graphviz label for a designator to an output stream.
inline void graphvizLabelTags(const PortalTagDesignatorVariant& d,
                              std::ostream& os) {
  std::visit(
      [&]<typename T>(const T& des) {
        if constexpr (std::is_same_v<T, std::monostate>) {
          os << "<br/><i>NullPortalDesignator</i>";
        } else {
          des.graphvizLabel(os);
        }
      },
      d);
}

}  // namespace Acts::Experimental::detail
