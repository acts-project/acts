// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/MaterialMapper.hpp"

#include "Acts/Utilities/Enumerate.hpp"

Acts::MaterialMapper::MaterialMapper(const Config& cfg,
                                     std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}

std::unique_ptr<Acts::MaterialMapper::State> Acts::MaterialMapper::createState()
    const {
  // Create the state
  auto state = std::make_unique<State>();
  // Create the surface material accumulater state
  state->surfaceMaterialAccumulaterState =
      m_cfg.surfaceMaterialAccumulater->createState();

  return state;
}

std::pair<Acts::RecordedMaterialTrack, Acts::RecordedMaterialTrack>
Acts::MaterialMapper::mapMaterial(State& state, const GeometryContext& gctx,
                                  const MagneticFieldContext& mctx,
                                  const RecordedMaterialTrack& rmTrack,
                                  const Options& options) const {
  // The recorded material track
  const auto& [starDir, recordedMaterial] = rmTrack;
  const auto& [position, direction] = starDir;
  auto [surfaceAssignments, volumeAssignments] =
      m_cfg.assignmentFinder->assignmentCandidates(gctx, mctx, position,
                                                   direction);

  // The mapped and unmapped material
  RecordedMaterialTrack mappedMaterial = {starDir, {}};
  RecordedMaterialTrack unmappedMaterial = {starDir, {}};
  // Assign the surface interactions
  auto [assigned, unassigned, emptyBinSurfaces] =
      MaterialInteractionAssignment::assign(
          gctx, recordedMaterial.materialInteractions, surfaceAssignments,
          options.assignmentOptions);

  // Record the assigned ones - as mapped ones
  mappedMaterial.second.materialInteractions.insert(
      mappedMaterial.second.materialInteractions.end(), assigned.begin(),
      assigned.end());

  // Record the unassigned ones - as unmapped ones
  unmappedMaterial.second.materialInteractions.insert(
      unmappedMaterial.second.materialInteractions.end(), unassigned.begin(),
      unassigned.end());

  // The material interactions
  m_cfg.surfaceMaterialAccumulater->accumulate(
      *state.surfaceMaterialAccumulaterState.get(), assigned, emptyBinSurfaces);

  return {mappedMaterial, unmappedMaterial};
}

Acts::MaterialMapper::DetectorMaterialMaps Acts::MaterialMapper::finalizeMaps(
    const State& state) const {
  // The final maps
  DetectorMaterialMaps detectorMaterialMaps;
  // The surface maps
  detectorMaterialMaps.first =
      m_cfg.surfaceMaterialAccumulater->finalizeMaterial(
          *state.surfaceMaterialAccumulaterState.get());

  return detectorMaterialMaps;
}