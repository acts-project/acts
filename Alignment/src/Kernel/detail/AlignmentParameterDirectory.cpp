// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsAlignment/Kernel/detail/AlignmentParameterDirectory.hpp"

#include <array>
#include <unordered_set>
#include <utility>

namespace ActsAlignment::detail {
namespace {

// AlignmentMask::All includes unused bits; enumerate only the six components.
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

using StructureSet = std::unordered_set<const AlignableStructure*>;

// Collect the hierarchy and reject repeated structures or member placements.
// The visited set also stops cycles before revisiting their children.
Acts::Result<StructureSet, AlignmentSelectionError> collectStructures(
    std::span<const std::shared_ptr<AlignableStructure>> roots) {
  StructureSet registered;
  std::unordered_set<const Acts::SurfacePlacementBase*> placements;
  std::vector<const AlignableStructure*> pending;
  for (const auto& root : roots) {
    pending.push_back(root.get());
  }
  while (!pending.empty()) {
    const auto* structure = pending.back();
    pending.pop_back();
    if (structure == nullptr) {
      return AlignmentSelectionError{AlignmentError::HierarchyValidationFailure,
                                     "Null hierarchy root or child"};
    }
    if (!registered.insert(structure).second) {
      return AlignmentSelectionError{AlignmentError::HierarchyValidationFailure,
                                     "Cycle or repeated structure ownership"};
    }
    for (const auto& surface : structure->surfacePointers()) {
      if (!surface) {
        return AlignmentSelectionError{
            AlignmentError::HierarchyValidationFailure,
            "Null structure member"};
      }
      const auto* placement = surface->surfacePlacement();
      if (placement != nullptr && !placements.insert(placement).second) {
        return AlignmentSelectionError{
            AlignmentError::HierarchyValidationFailure,
            "Repeated registered placement"};
      }
    }
    for (const auto& child : structure->children()) {
      pending.push_back(child.get());
    }
  }
  return registered;
}

}  // namespace

std::size_t AlignmentParameterDirectory::addAlignable(
    std::variant<const Acts::SurfacePlacementBase*, const AlignableStructure*>
        alignable,
    AlignmentMask mask) {
  AlignableUnit unit{alignable, {}};
  for (const auto& [dof, bit] : namedDofs) {
    if (ACTS_CHECK_BIT(mask, bit)) {
      unit.parameterIndices[dof] = numParameters++;
    }
  }
  const auto position = alignables.size();
  alignables.push_back(unit);
  return position;
}

bool AlignmentParameterDirectory::registerPlacement(
    const Acts::SurfacePlacementBase* placement, std::size_t alignableIndex) {
  return placementToAlignable.emplace(placement, alignableIndex).second;
}

Acts::Result<AlignmentParameterDirectory, AlignmentSelectionError>
makeAlignmentParameterDirectory(
    std::span<const std::shared_ptr<AlignableStructure>> structures,
    std::span<const std::shared_ptr<AlignableStructure>> selectedStructures,
    std::span<Acts::SurfacePlacementBase* const> modules,
    AlignmentMask moduleMask) {
  auto registered = collectStructures(structures);
  if (!registered.ok()) {
    return registered.error();
  }

  AlignmentParameterDirectory directory;
  for (const auto* module : modules) {
    if (module == nullptr) {
      return AlignmentSelectionError{AlignmentError::InvalidAlignmentSelection,
                                     "Null standalone module"};
    }
    const auto position = directory.addAlignable(module, moduleMask);
    if (!directory.registerPlacement(module, position)) {
      return AlignmentSelectionError{AlignmentError::InvalidAlignmentSelection,
                                     "Duplicate standalone module"};
    }
  }

  // A valid hierarchy can still be accompanied by a duplicate selection.
  StructureSet selected;
  for (const auto& structure : selectedStructures) {
    if (!structure || !registered->contains(structure.get())) {
      return AlignmentSelectionError{
          AlignmentError::InvalidAlignmentSelection,
          "Selected structure is null or unregistered"};
    }
    if (!selected.insert(structure.get()).second) {
      return AlignmentSelectionError{AlignmentError::InvalidAlignmentSelection,
                                     "Duplicate structure selection"};
    }
    const auto mask = structure->alignmentMask();
    if (!hasParameters(mask)) {
      continue;
    }
    if (!structure->children().empty() ||
        structure->surfacePointers().empty()) {
      return AlignmentSelectionError{
          AlignmentError::InvalidAlignmentSelection,
          "Active structure must be a nonempty leaf"};
    }
    const auto position = directory.addAlignable(structure.get(), mask);
    // Members share one entry. Registering them also detects overlap with
    // the standalone modules already in the lookup.
    for (const auto& surface : structure->surfacePointers()) {
      const auto* placement = surface->surfacePlacement();
      if (placement == nullptr) {
        return AlignmentSelectionError{
            AlignmentError::InvalidAlignmentSelection,
            "Selected member has no placement"};
      }
      if (!directory.registerPlacement(placement, position)) {
        return AlignmentSelectionError{
            AlignmentError::InvalidAlignmentSelection,
            "Overlapping alignment selections"};
      }
    }
  }
  return directory;
}

}  // namespace ActsAlignment::detail
