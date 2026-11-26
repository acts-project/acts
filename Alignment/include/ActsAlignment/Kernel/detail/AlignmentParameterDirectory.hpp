// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsAlignment/Geometry/AlignableStructure.hpp"
#include "ActsAlignment/Kernel/AlignmentError.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace ActsAlignment::detail {

/// The reason an alignment selection was rejected.
struct AlignmentSelectionError {
  AlignmentError code;
  std::string message;
};

/// Parameter indices for one selected module or structure.
struct AlignableUnit {
  /// Non-owning pointer to a selected module or structure.
  std::variant<const Acts::SurfacePlacementBase*, const AlignableStructure*>
      alignable;
  /// Position in the correction vector for each Acts::AlignmentIndices
  /// component. A fixed component has no index (std::nullopt).
  std::array<std::optional<std::size_t>, Acts::eAlignmentSize> parameterIndices;
};

/// Maps selected placements to alignables and their parameter indices.
struct AlignmentParameterDirectory {
  /// Number a unit's enabled parameters and append it to alignables.
  /// @param unit Module or structure to add; existing parameter indices are replaced
  /// @param mask Enabled translations and rotations
  /// @return Position of the new entry in the alignables vector
  std::size_t addAlignable(AlignableUnit unit, AlignmentMask mask);

  /// Associate a module placement with an alignable entry.
  /// @param placement Module placement to register
  /// @param alignableIndex Entry position returned by addAlignable()
  /// @return True if inserted; false if already registered, leaving its entry unchanged
  bool registerPlacement(const Acts::SurfacePlacementBase& placement,
                         std::size_t alignableIndex);

  /// Standalone modules in input order, followed by active selected structures.
  std::vector<AlignableUnit> alignables;
  /// Maps each module placement to its position in the alignables vector.
  /// A standalone module maps to its own entry. Every module belonging to a
  /// selected structure maps to that structure's shared entry.
  std::unordered_map<const Acts::SurfacePlacementBase*, std::size_t>
      placementToAlignable;
  /// Number of enabled parameters across all alignables.
  std::size_t numParameters = 0;
};

/// Number the enabled translations and rotations of selected alignables.
///
/// Active structures must have members and no children, and every member must
/// have a surface placement. Structures with mask None are omitted. Each
/// selected placement must belong to only one alignable, including standalone
/// modules with mask None.
///
/// Parameters are numbered from zero in input order: modules, then structures.
/// Within each alignable, Center0..2 precede Rotation0..2; fixed components are
/// skipped. An empty directory is valid. Rebuild after changing selection,
/// membership, or masks. Do not change the hierarchy while building the
/// directory. All selected structures and placements must outlive the
/// directory.
/// @param structures Top-level structures in the registered hierarchy
/// @param selectedStructures Registered structures to align
/// @param modules Modules selected for individual alignment, excluding members
/// of active selected structures
/// @param moduleMask Enabled translations and rotations for standalone modules
/// @return Parameter directory, or the reason the hierarchy or selection is invalid
Acts::Result<AlignmentParameterDirectory, AlignmentSelectionError>
makeAlignmentParameterDirectory(
    std::span<const std::shared_ptr<AlignableStructure>> structures,
    std::span<const std::shared_ptr<AlignableStructure>> selectedStructures,
    std::span<Acts::SurfacePlacementBase* const> modules,
    AlignmentMask moduleMask);

}  // namespace ActsAlignment::detail
