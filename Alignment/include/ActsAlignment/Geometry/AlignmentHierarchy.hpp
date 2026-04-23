// This file is part of the ACTS project.
//
// Copyright (C) 2026 CERN for the benefit of the ACTS project
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
#include <unordered_set>
#include <vector>

namespace ActsAlignment {

/// @brief Registry over a set of AlignableStructures with validation helpers.
///
/// AlignmentHierarchy builds a detector-element → owning-structure lookup at
/// construction and exposes it for later steps of the alignment pipeline
/// (validation, structure-level derivative accumulation). It does not own the
/// structures; the caller retains ownership.
class AlignmentHierarchy {
 public:
  /// @brief Outcome of validate()
  struct ValidationResult {
    /// Detector elements that appear in more than one structure. Each offender
    /// is listed once.
    std::vector<const Acts::SurfacePlacementBase*> overlapping;

    /// @brief Whether the hierarchy passes validation
    bool ok() const { return overlapping.empty(); }
  };

  /// @brief Construct from a flat list of structures
  /// @param structures Shared pointers to the alignable structures
  explicit AlignmentHierarchy(
      const std::vector<std::shared_ptr<AlignableStructure>>& structures)
      : m_structures(structures) {
    std::unordered_map<const Acts::SurfacePlacementBase*, unsigned int> counts;
    for (const auto& structure : m_structures) {
      if (structure == nullptr) {
        continue;
      }
      for (const auto& surface : structure->surfaces()) {
        const auto* detElement = surface->surfacePlacement();
        if (detElement == nullptr) {
          continue;
        }
        m_detElementToStructure.emplace(detElement, structure.get());
        counts[detElement]++;
      }
    }
    for (const auto& [detElement, count] : counts) {
      if (count > 1) {
        m_overlapping.push_back(detElement);
      }
    }
  }

  /// @brief Find the structure that owns a given detector element
  /// @param detElement The detector element to look up
  /// @return The owning structure, or nullptr if the element is a standalone
  ///         floating module
  AlignableStructure* structureFor(
      const Acts::SurfacePlacementBase* detElement) const {
    auto it = m_detElementToStructure.find(detElement);
    return it == m_detElementToStructure.end() ? nullptr : it->second;
  }

  /// @brief Check that no detector element is assigned to more than one
  ///        structure. Mixed mode is allowed: a detector element may be absent
  ///        from every structure and float as a standalone module.
  /// @return The validation result
  ValidationResult validate() const {
    return ValidationResult{m_overlapping};
  }

  /// @brief A DoF conflict between a structure and one of its child modules.
  struct MaskConflict {
    /// The offending detector element
    const Acts::SurfacePlacementBase* detElement;
    /// The parent structure whose mask overlaps
    AlignableStructure* structure;
    /// The bits that are set in both masks
    AlignmentMask conflictingBits;
  };

  /// @brief Detect DoF bits that are simultaneously floating at structure and
  ///        module level for the same detector element.
  ///
  /// A structure-level DoF and a module-level DoF on the same physical
  /// element are correlated: if both are floating, the Hessian picks up a
  /// zero eigenvalue along that direction. This check returns every such
  /// overlap so the caller can reject the configuration.
  /// @param moduleMask The DoFs currently floating at module level
  ///                   (i.e. the iteration's global alignment mask)
  /// @param floatingModules Detector elements listed in @c alignedDetElements
  /// @return One entry per conflicting (detElement, structure) pair
  std::vector<MaskConflict> detectMaskConflicts(
      AlignmentMask moduleMask,
      const std::vector<Acts::SurfacePlacementBase*>& floatingModules) const {
    std::unordered_set<const Acts::SurfacePlacementBase*> floatingSet(
        floatingModules.begin(), floatingModules.end());
    std::vector<MaskConflict> conflicts;
    for (const auto& [detElement, structure] : m_detElementToStructure) {
      if (!floatingSet.contains(detElement)) {
        continue;
      }
      AlignmentMask overlap = structure->alignmentMask() & moduleMask;
      if (overlap != AlignmentMask::None) {
        conflicts.push_back({detElement, structure, overlap});
      }
    }
    return conflicts;
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
  std::unordered_map<const Acts::SurfacePlacementBase*, AlignableStructure*>
      m_detElementToStructure;

  /// Detector elements that appear in more than one structure
  std::vector<const Acts::SurfacePlacementBase*> m_overlapping;
};

}  // namespace ActsAlignment
