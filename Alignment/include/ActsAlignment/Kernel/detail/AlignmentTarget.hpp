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
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
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

/// One optional compact column per named alignment DoF; nullopt means fixed.
using AlignmentParameterIndices =
    std::array<std::optional<std::size_t>, Acts::eAlignmentSize>;

/// Reverse mapping from a compact column to its target and named DoF.
struct AlignmentParameter {
  std::size_t target;
  Acts::AlignmentIndices dof;
};

/// Iteration-specific compact parameter layout, built from a prepared
/// selection.
struct AlignmentParameterLayout {
  /// Indexed by target in PreparedAlignmentTargets, including fixed modules.
  std::vector<AlignmentParameterIndices> targetColumns;
  /// Compact columns in target order, then Center0..2, Rotation0..2 order.
  std::vector<AlignmentParameter> parameters;
  /// Only placements with active columns appear here.
  std::unordered_map<const Acts::SurfacePlacementBase*, std::size_t>
      placementTargets;
};

/// Measurement rows supplied to the track layout builder.
struct AlignmentMeasurement {
  /// Null for rows without an associated placement (e.g. track-model rows).
  const Acts::SurfacePlacementBase* placement = nullptr;
  std::size_t dimension = 0;
};

/// One target present on a track, with indices into its local parameter vector.
struct TrackAlignmentTarget {
  std::size_t globalTarget = 0;
  AlignmentParameterIndices columns;
};

/// Measurement metadata is independent of the number of active parameters.
struct TrackAlignmentMeasurement {
  /// Index into TrackAlignmentLayout::targets, or nullopt for fixed/target-less
  /// rows.
  std::optional<std::size_t> target;
  std::size_t rowOffset = 0;
  std::size_t dimension = 0;
};

/// Mapping only: does not compute residuals, derivatives, or a track-model
/// gradient.
struct TrackAlignmentLayout {
  /// Unique participating targets in global target order.
  std::vector<TrackAlignmentTarget> targets;
  /// Global column for each local column, for normal-equation accumulation.
  std::vector<std::size_t> localToGlobal;
  /// One entry for every input measurement, including repeated member hits.
  std::vector<TrackAlignmentMeasurement> measurements;
  std::size_t measurementDim = 0;
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
/// @param activeStructures Explicitly selected registered structures. Structures
/// with active DoFs must be leaves; None selections stay metadata and may have
/// children.
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

/// Build a compact layout with module masks resolved for this exact iteration.
/// Structure masks come from the snapshot and are independent of the schedule.
/// @param selection Snapshot returned by prepareAlignmentTargets
/// @param iteration Current iteration number
/// @param iterationState Module masks; a missing entry means All
/// @return Compact layout, or NoActiveAlignmentParameters if it has no columns
Acts::Result<AlignmentParameterLayout, AlignmentTargetError>
makeAlignmentParameterLayout(
    const PreparedAlignmentTargets& selection, unsigned int iteration,
    const std::map<unsigned int, AlignmentMask>& iterationState = {});

/// Map a track's measurement rows to shared columns. Scratch storage scales
/// with the track's measurements, independently of the global target count.
/// @param layout Global layout returned by makeAlignmentParameterLayout
/// @param measurements Rows in track order; may include fixed/target-less rows
/// @return Track-local mapping (possibly with no active targets)
TrackAlignmentLayout makeTrackAlignmentLayout(
    const AlignmentParameterLayout& layout,
    std::span<const AlignmentMeasurement> measurements);

}  // namespace ActsAlignment::detail
