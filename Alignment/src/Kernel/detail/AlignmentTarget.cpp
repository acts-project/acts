// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsAlignment/Kernel/detail/AlignmentTarget.hpp"

#include "Acts/Surfaces/SurfacePlacementBase.hpp"

#include <array>
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

}  // namespace ActsAlignment::detail
