// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/MaterialMapper.hpp"

#include "Acts/Utilities/detail/WorkerPool.hpp"

Acts::MaterialMapper::MaterialMapper(const Config& cfg,
                                     std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {
  if (m_cfg.assignmentFinder == nullptr) {
    throw std::invalid_argument("The assignment finder is not set");
  }
  if (m_cfg.surfaceMaterialAccumulator == nullptr) {
    throw std::invalid_argument("The surface material accumulator is not set");
  }
  if (m_cfg.numThreads > 1) {
    m_pool = std::make_shared<detail::WorkerPool>(m_cfg.numThreads);
  }
}

std::unique_ptr<Acts::MaterialMapper::State> Acts::MaterialMapper::createState(
    const GeometryContext& gctx) const {
  // Create the state
  auto state = std::make_unique<State>();
  // Create the surface material accumulator state
  state->surfaceMaterialAccumulatorState =
      m_cfg.surfaceMaterialAccumulator->createState(gctx);
  // Return the state object
  return state;
}

Acts::MaterialMapper::TrackAssignment Acts::MaterialMapper::assignTrack(
    const GeometryContext& gctx, const MagneticFieldContext& mctx,
    const RecordedMaterialTrack& rmTrack, const Options& options) const {
  // The recorded material track
  const auto& [starDir, recordedMaterial] = rmTrack;
  const auto& [position, direction] = starDir;
  auto [surfaceAssignments, volumeAssignments] =
      m_cfg.assignmentFinder->assignmentCandidates(gctx, mctx, position,
                                                   direction);

  // The mapped and unmapped material
  TrackAssignment result{{starDir, {}}, {starDir, {}}, {}};
  // Assign the surface interactions
  auto [assigned, unassigned, emptyBinSurfaces] =
      MaterialInteractionAssignment::assign(
          gctx, recordedMaterial.materialInteractions, surfaceAssignments,
          options.assignmentOptions);

  // Record the assigned ones - as mapped ones
  result.mapped.second.materialInteractions = std::move(assigned);
  // Record the unassigned ones - as unmapped ones
  result.unmapped.second.materialInteractions = std::move(unassigned);
  result.surfacesWithoutAssignment = std::move(emptyBinSurfaces);

  // The function to calculate the total material before returning
  auto calculateTotalMaterial = [](RecordedMaterialTrack& rTrack) -> void {
    for (const auto& mi : rTrack.second.materialInteractions) {
      rTrack.second.materialInX0 += mi.materialSlab.thicknessInX0();
      rTrack.second.materialInL0 += mi.materialSlab.thicknessInL0();
    }
  };
  // Fill the totals to the material tracks (useful for debugging)
  calculateTotalMaterial(result.mapped);
  calculateTotalMaterial(result.unmapped);
  return result;
}

std::pair<Acts::RecordedMaterialTrack, Acts::RecordedMaterialTrack>
Acts::MaterialMapper::mapMaterial(State& state, const GeometryContext& gctx,
                                  const MagneticFieldContext& mctx,
                                  const RecordedMaterialTrack& rmTrack,
                                  const Options& options) const {
  // A batch of one track
  const RecordedMaterialTrack* track = &rmTrack;
  return std::move(
      mapMaterial(state, gctx, mctx, std::span(&track, 1), options).front());
}

std::vector<std::pair<Acts::RecordedMaterialTrack, Acts::RecordedMaterialTrack>>
Acts::MaterialMapper::mapMaterial(
    State& state, const GeometryContext& gctx, const MagneticFieldContext& mctx,
    std::span<const RecordedMaterialTrack* const> rmTracks,
    const Options& options) const {
  // The assignment of the tracks is independent: run it concurrently
  std::vector<TrackAssignment> assignments(rmTracks.size());
  auto assign = [&](std::size_t i) {
    assignments[i] = assignTrack(gctx, mctx, *rmTracks[i], options);
  };
  if (m_pool != nullptr) {
    m_pool->parallelFor(rmTracks.size(), assign);
  } else {
    for (std::size_t i = 0; i < rmTracks.size(); ++i) {
      assign(i);
    }
  }

  // The accumulation follows the track order, the accumulator decides
  // which parts of it can run concurrently
  std::vector<ISurfaceMaterialAccumulator::TrackAccumulation> accumulations;
  accumulations.reserve(assignments.size());
  for (const auto& assignment : assignments) {
    accumulations.push_back({&assignment.mapped.second.materialInteractions,
                             &assignment.surfacesWithoutAssignment});
  }
  m_cfg.surfaceMaterialAccumulator->accumulateTracks(
      *state.surfaceMaterialAccumulatorState, gctx, accumulations,
      m_pool.get());

  std::vector<std::pair<RecordedMaterialTrack, RecordedMaterialTrack>> result;
  result.reserve(assignments.size());
  for (auto& assignment : assignments) {
    result.emplace_back(std::move(assignment.mapped),
                        std::move(assignment.unmapped));
  }
  return result;
}

Acts::TrackingGeometryMaterial Acts::MaterialMapper::finalizeMaps(
    const State& state, const GeometryContext& gctx) const {
  return m_cfg.surfaceMaterialAccumulator->finalizeMaps(
      *state.surfaceMaterialAccumulatorState, gctx);
}
