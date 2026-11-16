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

#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace ActsAlignment::detail {

/// Preparation failure with diagnostic context, without logging an input error.
struct AlignmentTargetError {
  AlignmentError code;
  /// Target/phase and original provider message, when available.
  std::string message;
};

/// One placement and its initial contextual transform. The placement is
/// borrowed.
struct AlignmentTargetMember {
  const Acts::SurfacePlacementBase* placement = nullptr;
  Acts::Transform3 transform = Acts::Transform3::Identity();
};

/// A numerical target, independent of track measurements and geometry IDs.
struct AlignmentTarget {
  enum class Kind { Module, Structure };
  Kind kind = Kind::Module;
  /// Retains structure ownership; null for a standalone module.
  std::shared_ptr<const AlignableStructure> structure;
  /// Snapshotted direct membership; callers retain all placement lifetimes.
  std::vector<AlignmentTargetMember> members;
  /// Initial frame, present only for a structure. Later updates use working
  /// state.
  std::optional<Acts::Transform3> frame;
  /// Captured structure mask. Module masks are supplied per iteration.
  AlignmentMask mask = AlignmentMask::All;
};

/// Per-call snapshot. Building layouts does not query providers again.
/// Do not mutate declarations or geometry concurrently with preparation.
struct PreparedAlignmentTargets {
  /// Retain registered structures, including inactive metadata.
  std::vector<std::shared_ptr<AlignableStructure>> structures;
  /// Modules in input order followed by non-None structures in selection order.
  std::vector<AlignmentTarget> targets;
};

/// Read-only declaration of whether an application can update a placement.
/// Must not publish transforms. This does not guarantee transactional updates;
/// the eventual numerical solver must also validate its publication interface.
using AlignmentUpdateSupport =
    std::function<bool(const Acts::SurfacePlacementBase&)>;

/// Validate the graph, registered placement uniqueness, and explicit selection.
/// Safe before constructing AlignmentHierarchy, including for cyclic input.
/// None selections require no numerical frame and may have
/// children/constraints.
/// @param structures Registered roots
/// @param activeStructures Explicit selections, which must be registered
/// @return Whether any selected structure has active DoFs, or a diagnostic
Acts::Result<bool, AlignmentTargetError> validateAlignmentSelection(
    std::span<const std::shared_ptr<AlignableStructure>> structures,
    std::span<const std::shared_ptr<AlignableStructure>> activeStructures);

/// Prepare targets without fitting or changing geometry. Resolve frames once.
/// Validates proper rigid transforms with absolute tolerance 1e-10 on the
/// orthogonality, determinant and affine-row tests. Translation units are ACTS
/// length units; all structure DoFs use its local frame.
/// Recoverable provider/capability exceptions become diagnostic failures;
/// std::bad_alloc and non-standard exceptions propagate.
/// @param gctx Initial fit geometry context, also passed to member placements
/// @param structures Registered roots, retained by the returned selection
/// @param activeStructures Explicit leaf selection; None entries stay metadata
/// @param modules Standalone placements; callers retain their lifetimes
/// @param supportsUpdate Read-only per-placement update capability; required for targets
/// @return A call-owned snapshot or a diagnostic (an empty snapshot is allowed)
Acts::Result<PreparedAlignmentTargets, AlignmentTargetError>
prepareAlignmentTargets(
    const Acts::GeometryContext& gctx,
    std::span<const std::shared_ptr<AlignableStructure>> structures,
    std::span<const std::shared_ptr<AlignableStructure>> activeStructures,
    std::span<Acts::SurfacePlacementBase* const> modules,
    const AlignmentUpdateSupport& supportsUpdate);

}  // namespace ActsAlignment::detail
