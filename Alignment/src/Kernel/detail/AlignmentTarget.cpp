// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsAlignment/Kernel/detail/AlignmentTarget.hpp"

#include "Acts/Surfaces/SurfacePlacementBase.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <new>
#include <sstream>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ActsAlignment::detail {
namespace {

constexpr std::array<std::pair<Acts::AlignmentIndices, AlignmentMask>,
                     Acts::eAlignmentSize>
    namedDofs{{{Acts::eAlignmentCenter0, AlignmentMask::Center0},
               {Acts::eAlignmentCenter1, AlignmentMask::Center1},
               {Acts::eAlignmentCenter2, AlignmentMask::Center2},
               {Acts::eAlignmentRotation0, AlignmentMask::Rotation0},
               {Acts::eAlignmentRotation1, AlignmentMask::Rotation1},
               {Acts::eAlignmentRotation2, AlignmentMask::Rotation2}}};

bool hasParameters(AlignmentMask mask) {
  for (const auto& [dof, bit] : namedDofs) {
    if (ACTS_CHECK_BIT(mask, bit)) {
      return true;
    }
  }
  return false;
}

std::string describe(const AlignableStructure& structure) {
  std::ostringstream os;
  os << "structure " << structure.geometryId() << " (" << &structure << ")";
  return os.str();
}

bool isRigid(const Acts::Transform3& transform) {
  constexpr double tolerance = 1e-10;
  return transform.matrix().allFinite() &&
         (transform.matrix().row(3) - Acts::Vector4(0, 0, 0, 1).transpose())
                 .cwiseAbs()
                 .maxCoeff() <= tolerance &&
         (transform.linear().transpose() * transform.linear() -
          Acts::RotationMatrix3::Identity())
                 .cwiseAbs()
                 .maxCoeff() <= tolerance &&
         std::abs(transform.linear().determinant() - 1.) <= tolerance;
}

using StructureSet = std::unordered_set<const AlignableStructure*>;

Acts::Result<StructureSet, AlignmentTargetError> inspectHierarchy(
    std::span<const std::shared_ptr<AlignableStructure>> structures) {
  StructureSet registered;
  std::unordered_set<const Acts::SurfacePlacementBase*> placements;
  std::vector<const AlignableStructure*> pending;
  for (auto it = structures.rbegin(); it != structures.rend(); ++it) {
    pending.push_back(it->get());
  }
  while (!pending.empty()) {
    const auto* structure = pending.back();
    pending.pop_back();
    if (structure == nullptr) {
      return AlignmentTargetError{AlignmentError::HierarchyValidationFailure,
                                  "Null hierarchy root or child"};
    }
    if (!registered.insert(structure).second) {
      return AlignmentTargetError{
          AlignmentError::HierarchyValidationFailure,
          "Cycle or repeated ownership of " + describe(*structure)};
    }
    for (const auto& surface : structure->surfacePointers()) {
      if (!surface) {
        return AlignmentTargetError{AlignmentError::HierarchyValidationFailure,
                                    "Null member of " + describe(*structure)};
      }
      const auto* placement = surface->surfacePlacement();
      if (placement != nullptr && !placements.insert(placement).second) {
        return AlignmentTargetError{
            AlignmentError::HierarchyValidationFailure,
            "Repeated registered placement in " + describe(*structure)};
      }
    }
    const auto children = structure->children();
    for (auto it = children.rbegin(); it != children.rend(); ++it) {
      pending.push_back(it->get());
    }
  }
  return registered;
}

}  // namespace

Acts::Result<bool, AlignmentTargetError> validateAlignmentSelection(
    std::span<const std::shared_ptr<AlignableStructure>> structures,
    std::span<const std::shared_ptr<AlignableStructure>> activeStructures) {
  auto registered = inspectHierarchy(structures);
  if (!registered.ok()) {
    return registered.error();
  }
  StructureSet selected;
  bool hasActiveStructures = false;
  for (std::size_t i = 0; i < activeStructures.size(); ++i) {
    const auto& structure = activeStructures[i];
    const std::string label = "Selected structure " + std::to_string(i);
    if (!structure || !registered->contains(structure.get())) {
      return AlignmentTargetError{AlignmentError::InvalidAlignmentSelection,
                                  label + " is null or unregistered"};
    }
    if (!selected.insert(structure.get()).second) {
      return AlignmentTargetError{
          AlignmentError::InvalidAlignmentSelection,
          "Duplicate selection of " + describe(*structure)};
    }
    if (!hasParameters(structure->alignmentMask())) {
      continue;
    }
    if (!structure->children().empty()) {
      return AlignmentTargetError{AlignmentError::InvalidAlignmentSelection,
                                  "Active non-leaf " + describe(*structure)};
    }
    hasActiveStructures = true;
  }
  return hasActiveStructures;
}

Acts::Result<PreparedAlignmentTargets, AlignmentTargetError>
prepareAlignmentTargets(
    const Acts::GeometryContext& gctx,
    std::span<const std::shared_ptr<AlignableStructure>> structures,
    std::span<const std::shared_ptr<AlignableStructure>> activeStructures,
    std::span<Acts::SurfacePlacementBase* const> modules,
    const AlignmentUpdateSupport& supportsUpdate) {
  auto validation = validateAlignmentSelection(structures, activeStructures);
  if (!validation.ok()) {
    return validation.error();
  }
  PreparedAlignmentTargets prepared;
  prepared.structures.assign(structures.begin(), structures.end());
  for (const auto* module : modules) {
    if (module == nullptr) {
      return AlignmentTargetError{AlignmentError::InvalidAlignmentSelection,
                                  "Null standalone module"};
    }
    AlignmentTarget target;
    target.members.push_back({module, Acts::Transform3::Identity()});
    prepared.targets.push_back(std::move(target));
  }
  for (const auto& structure : activeStructures) {
    const auto mask = structure->alignmentMask();
    if (!hasParameters(mask)) {
      continue;
    }
    if (!structure->constraints().empty()) {
      return AlignmentTargetError{
          AlignmentError::UnsupportedAlignmentConstraints,
          "Constraints on " + describe(*structure)};
    }
    if (structure->surfacePointers().empty()) {
      return AlignmentTargetError{AlignmentError::InvalidAlignmentSelection,
                                  "Empty active " + describe(*structure)};
    }
    AlignmentTarget target;
    target.kind = AlignmentTarget::Kind::Structure;
    target.structure = structure;
    target.mask = mask;
    for (const auto& surface : structure->surfacePointers()) {
      const auto* placement = surface->surfacePlacement();
      if (placement == nullptr || &placement->surface() != surface.get()) {
        return AlignmentTargetError{
            AlignmentError::InvalidAlignmentSelection,
            "Missing placement or invalid member backlink in " +
                describe(*structure)};
      }
      target.members.push_back({placement, Acts::Transform3::Identity()});
    }
    prepared.targets.push_back(std::move(target));
  }

  std::unordered_set<const Acts::SurfacePlacementBase*> claimed;
  for (std::size_t i = 0; i < prepared.targets.size(); ++i) {
    auto& target = prepared.targets[i];
    const std::string label = target.structure ? describe(*target.structure)
                                               : "module " + std::to_string(i);
    for (const auto& member : target.members) {
      if (!claimed.insert(member.placement).second) {
        return AlignmentTargetError{
            AlignmentError::InvalidAlignmentSelection,
            "Duplicate/overlapping placement at " + label};
      }
      if (member.placement->surface().surfacePlacement() != member.placement) {
        return AlignmentTargetError{AlignmentError::InvalidAlignmentSelection,
                                    "Invalid placement backlink at " + label};
      }
    }
  }

  // All structural checks precede provider/capability calls. These calls only
  // read the caller's state; nothing is published, even on partial failure.
  for (std::size_t i = 0; i < prepared.targets.size(); ++i) {
    auto& target = prepared.targets[i];
    const std::string label = target.structure ? describe(*target.structure)
                                               : "module " + std::to_string(i);
    AlignmentError code = AlignmentError::InvalidAlignmentSelection;
    std::string phase = "update capability";
    try {
      for (const auto& member : target.members) {
        if (!supportsUpdate || !supportsUpdate(*member.placement)) {
          return AlignmentTargetError{code, label + ": no update route"};
        }
      }
      code = AlignmentError::AlignmentFrameResolutionFailure;
      phase = "initial structure frame in fit geometry context";
      if (target.structure) {
        target.frame = target.structure->localToGlobalTransform(gctx);
        if (!isRigid(*target.frame)) {
          return AlignmentTargetError{code, label + ": invalid " + phase};
        }
      }
      for (std::size_t j = 0; j < target.members.size(); ++j) {
        phase = "initial member " + std::to_string(j) +
                " frame in fit geometry context";
        auto& member = target.members[j];
        member.transform = member.placement->localToGlobalTransform(gctx);
        if (!isRigid(member.transform)) {
          return AlignmentTargetError{code, label + ": invalid " + phase};
        }
      }
    } catch (const std::bad_alloc&) {
      throw;
    } catch (const std::exception& error) {
      return AlignmentTargetError{code,
                                  label + ": " + phase + ": " + error.what()};
    }
  }
  return prepared;
}

Acts::Result<AlignmentParameterLayout, AlignmentTargetError>
makeAlignmentParameterLayout(
    const PreparedAlignmentTargets& selection, unsigned int iteration,
    const std::map<unsigned int, AlignmentMask>& iterationState) {
  const auto it = iterationState.find(iteration);
  const auto moduleMask =
      it == iterationState.end() ? AlignmentMask::All : it->second;
  AlignmentParameterLayout layout;
  layout.targetColumns.resize(selection.targets.size());
  for (std::size_t i = 0; i < selection.targets.size(); ++i) {
    const auto& target = selection.targets[i];
    const auto mask =
        target.kind == AlignmentTarget::Kind::Module ? moduleMask : target.mask;
    for (const auto& [dof, bit] : namedDofs) {
      if (ACTS_CHECK_BIT(mask, bit)) {
        layout.targetColumns[i][dof] = layout.parameters.size();
        layout.parameters.push_back({i, dof});
      }
    }
    if (hasParameters(mask)) {
      for (const auto& member : target.members) {
        layout.placementTargets.emplace(member.placement, i);
      }
    }
  }
  if (layout.parameters.empty()) {
    return AlignmentTargetError{AlignmentError::NoActiveAlignmentParameters,
                                "Iteration " + std::to_string(iteration) +
                                    " has no active alignment parameters"};
  }
  return layout;
}

TrackAlignmentLayout makeTrackAlignmentLayout(
    const AlignmentParameterLayout& layout,
    std::span<const AlignmentMeasurement> measurements) {
  std::vector<std::size_t> participating;
  participating.reserve(measurements.size());
  for (const auto& measurement : measurements) {
    const auto it = layout.placementTargets.find(measurement.placement);
    if (it != layout.placementTargets.end()) {
      participating.push_back(it->second);
    }
  }
  std::sort(participating.begin(), participating.end());
  participating.erase(std::unique(participating.begin(), participating.end()),
                      participating.end());
  TrackAlignmentLayout track;
  track.targets.reserve(participating.size());
  track.measurements.reserve(measurements.size());
  for (const auto i : participating) {
    TrackAlignmentTarget target{i, {}};
    for (const auto& [dof, bit] : namedDofs) {
      if (const auto column = layout.targetColumns[i][dof]; column) {
        target.columns[dof] = track.localToGlobal.size();
        track.localToGlobal.push_back(*column);
      }
    }
    track.targets.push_back(target);
  }
  for (const auto& measurement : measurements) {
    std::optional<std::size_t> target;
    const auto it = layout.placementTargets.find(measurement.placement);
    if (it != layout.placementTargets.end()) {
      target = static_cast<std::size_t>(std::lower_bound(participating.begin(),
                                                         participating.end(),
                                                         it->second) -
                                        participating.begin());
    }
    track.measurements.push_back(
        {target, track.measurementDim, measurement.dimension});
    track.measurementDim += measurement.dimension;
  }
  return track;
}

}  // namespace ActsAlignment::detail
