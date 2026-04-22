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

#include <unordered_map>
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
  /// @param structures Non-owning pointers to the alignable structures
  explicit AlignmentHierarchy(
      const std::vector<AlignableStructure*>& structures)
      : m_structures(structures) {
    std::unordered_map<const Acts::SurfacePlacementBase*, unsigned int> counts;
    for (auto* structure : m_structures) {
      if (structure == nullptr) {
        continue;
      }
      for (const auto& surface : structure->surfaces()) {
        const auto* detElement = surface->surfacePlacement();
        if (detElement == nullptr) {
          continue;
        }
        m_detElementToStructure.emplace(detElement, structure);
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

  /// @brief Access the registered structures
  /// @return A vector of non-owning pointers to the structures
  const std::vector<AlignableStructure*>& structures() const {
    return m_structures;
  }

 private:
  /// The registered structures (non-owning)
  std::vector<AlignableStructure*> m_structures;

  /// Flat lookup from detector element to owning structure
  std::unordered_map<const Acts::SurfacePlacementBase*, AlignableStructure*>
      m_detElementToStructure;

  /// Detector elements that appear in more than one structure
  std::vector<const Acts::SurfacePlacementBase*> m_overlapping;
};

}  // namespace ActsAlignment
