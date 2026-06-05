// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsAlignment/Geometry/AlignableStructure.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

namespace ActsAlignment {

/// @brief Registry over a set of AlignableStructures with validation helpers.
///
/// AlignmentHierarchy builds a detector-element → owning-structure lookup at
/// construction and exposes it for later steps of the alignment pipeline
/// (validation, structure-level derivative accumulation). It does not own the
/// structures; the caller retains ownership.
///
/// @par Multi-level alignment
/// The hierarchy is rebuilt on every call to @c Alignment::align(), so
/// iterative multi-level alignment (strip system → endcap → disk → module,
/// one level active per call, cycling until convergence) is supported out of
/// the box. Each surface is mapped to its immediate owning structure; at most
/// one level of the hierarchy should have floating DoFs per @c align() call.
/// Simultaneous fitting of DoFs at multiple levels in a single call would
/// require extending the lookup to return the full ancestor chain and
/// accumulating chain-rule Jacobians at each level — this is not implemented.
class AlignmentHierarchy {
 public:
  /// @brief Outcome of validate()
  struct ValidationResult {
    /// Surfaces whose detector element appears in more than one structure.
    /// One representative surface is listed per conflicting element.
    std::vector<const Acts::Surface*> overlapping;

    /// @brief Whether the hierarchy passes validation
    bool ok() const { return overlapping.empty(); }
  };

  /// @brief Construct from a flat list of structures
  /// @param structures Shared pointers to the alignable structures
  explicit AlignmentHierarchy(
      const std::vector<std::shared_ptr<AlignableStructure>>& structures)
      : m_structures(structures) {
    std::unordered_map<const Acts::SurfacePlacementBase*, const Acts::Surface*>
        firstSurface;
    // Register a structure's direct surfaces, then recurse into children.
    // Each surface is registered under its immediate owning structure so that
    // child-level DoFs are not silently subsumed by the parent.
    auto registerStructure = [&](auto& self,
                                 AlignableStructure& structure) -> void {
      for (const Acts::Surface& surface : structure.surfaces()) {
        const auto* detElement = surface.surfacePlacement();
        if (detElement == nullptr) {
          continue;
        }
        const bool newElement =
            m_detElementToStructure.emplace(detElement, &structure).second;
        const auto surfIt = firstSurface.emplace(detElement, &surface).first;
        if (!newElement) {
          m_overlapping.push_back(surfIt->second);
        }
      }
      for (const auto& child : structure.children()) {
        if (child != nullptr) {
          self(self, *child);
        }
      }
    };
    for (const auto& structure : m_structures) {
      if (structure != nullptr) {
        registerStructure(registerStructure, *structure);
      }
    }
  }

  /// @brief Find the structure that owns a given detector element
  /// @param detElement The detector element to look up
  /// @return The owning structure, or nullptr if the element is a standalone
  ///         floating module
  const AlignableStructure* structureFor(
      const Acts::SurfacePlacementBase& detElement) const {
    auto it = m_detElementToStructure.find(&detElement);
    return it == m_detElementToStructure.end() ? nullptr : it->second;
  }

  /// @brief Find the structure that owns a given surface
  /// @param surface The surface to look up
  /// @return The owning structure, or nullptr if the surface's element is a
  ///         standalone floating module (or has no attached placement)
  const AlignableStructure* structureFor(const Acts::Surface& surface) const {
    const auto* detElement = surface.surfacePlacement();
    if (detElement == nullptr) {
      return nullptr;
    }
    return structureFor(*detElement);
  }

  /// @brief Check that no detector element is assigned to more than one
  ///        structure. Mixed mode is allowed: a detector element may be absent
  ///        from every structure and float as a standalone module.
  /// @return The validation result
  ValidationResult validate() const { return ValidationResult{m_overlapping}; }

  /// @brief A DoF conflict between a structure and one of its child modules.
  struct MaskConflict {
    /// The offending detector element
    const Acts::SurfacePlacementBase* detElement;
    /// The parent structure whose mask overlaps
    const AlignableStructure* structure;
    /// The bits that are set in both masks
    AlignmentMask conflictingBits;
  };

  /// @brief Detect DoF conflicts between structure-level and module-level
  ///        floating parameters.
  ///
  /// @todo Not yet implemented. A meaningful conflict check requires the
  ///       chain-rule Jacobian @f$\partial\alpha_\text{surface} /
  ///       \partial\alpha_\text{structure}@f$ to project structure-level DoFs
  ///       into the surface frame before comparing masks. Without it, a
  ///       bitwise mask comparison conflates DoF indices across different
  ///       coordinate frames and produces both false positives and false
  ///       negatives. Genuine redundancies are caught at solve time by a
  ///       rank-deficient Hessian.
  /// @return Always empty until the chain-rule Jacobian is available.
  std::vector<MaskConflict> detectMaskConflicts(
      AlignmentMask /*moduleMask*/,
      const std::vector<Acts::SurfacePlacementBase*>& /*floatingModules*/)
      const {
    return {};
  }

  /// @brief Access the registered structures
  /// @return A vector of shared pointers to the structures
  const std::vector<std::shared_ptr<AlignableStructure>>& structures() const {
    return m_structures;
  }

 private:
  /// The registered structures
  std::vector<std::shared_ptr<AlignableStructure>> m_structures;

  /// Flat lookup from detector element to owning structure
  std::unordered_map<const Acts::SurfacePlacementBase*,
                     const AlignableStructure*>
      m_detElementToStructure;

  /// One representative surface per detector element claimed by more than one
  /// structure
  std::vector<const Acts::Surface*> m_overlapping;
};

}  // namespace ActsAlignment
