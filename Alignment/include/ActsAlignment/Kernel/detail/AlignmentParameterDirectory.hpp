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
  /// Standalone modules in input order, followed by active selected structures.
  std::vector<AlignableUnit> alignables;
  /// Position in alignables for each selected placement.
  std::unordered_map<const Acts::SurfacePlacementBase*, std::size_t>
      placementToAlignable;
  /// Number of enabled parameters across all alignables.
  std::size_t numParameters = 0;
};

/// Number the enabled translations and rotations of selected alignables.
///
/// Active structures must have members and no children. Structures with mask
/// None are omitted. Each selected placement must belong to only one alignable,
/// including standalone modules with mask None.
///
/// Parameters are numbered from zero in input order: modules, then structures.
/// Within each alignable, Center0..2 precede Rotation0..2; fixed components are
/// skipped. An empty directory is valid. Rebuild after changing selection,
/// membership, or masks. Do not change the hierarchy while building the
/// directory. All selected structures and placements must outlive the
/// directory.
/// @param structures Top-level structures in the registered hierarchy
/// @param selectedStructures Registered structures to align
/// @param modules Standalone module placements
/// @param moduleMask Enabled translations and rotations for standalone modules
/// @return Parameter directory, or the reason the hierarchy or selection is invalid
Acts::Result<AlignmentParameterDirectory, AlignmentSelectionError>
makeAlignmentParameterDirectory(
    std::span<const std::shared_ptr<AlignableStructure>> structures,
    std::span<const std::shared_ptr<AlignableStructure>> selectedStructures,
    std::span<Acts::SurfacePlacementBase* const> modules,
    AlignmentMask moduleMask);

}  // namespace ActsAlignment::detail
