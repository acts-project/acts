// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/MaterialMapping/MaterialMapping.hpp"

#include "Acts/Material/AccumulatedMaterialSlab.hpp"
#include "Acts/Material/AccumulatedSurfaceMaterial.hpp"
#include "ActsExamples/MaterialMapping/IMaterialWriter.hpp"

#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace ActsExamples {

MaterialMapping::MaterialMapping(const MaterialMapping::Config& cfg,
                                 std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("MaterialMapping", std::move(logger)), m_cfg(cfg) {
  // Prepare the I/O collections
  m_inputMaterialTracks.initialize(m_cfg.inputMaterialTracks);
  m_outputMappedMaterialTracks.initialize(m_cfg.mappedMaterialTracks);
  m_outputUnmappedMaterialTracks.initialize(m_cfg.unmappedMaterialTracks);

  ACTS_LOG_WITH_LOGGER(this->logger(), Acts::Logging::INFO,
                       "This algorithm requires inter-event information, "
                           << "run in single-threaded mode! Use the "
                           << "numThreads option of the material mapper to "
                           << "map the tracks of each event concurrently.");

  if (m_cfg.materialMapper == nullptr) {
    throw std::invalid_argument("Missing material mapper");
  }
  // Create the state object
  m_mappingState = m_cfg.materialMapper->createState(m_cfg.geoContext);
}

ProcessCode MaterialMapping::finalize() {
  m_material =
      m_cfg.materialMapper->finalizeMaps(*m_mappingState, m_cfg.geoContext);
  for (auto& writer : m_cfg.materialWriters) {
    writer->writeMaterial(*m_material);
  }
  return ProcessCode::SUCCESS;
}

const Acts::TrackingGeometryMaterial& MaterialMapping::material() const {
  if (!m_material) {
    throw std::logic_error("Material maps are available after finalization");
  }
  return *m_material;
}

ProcessCode MaterialMapping::execute(const AlgorithmContext& context) const {
  // Take the collection from the EventStore: input collection
  const auto& mtrackCollection = m_inputMaterialTracks(context);

  // Write the output collections to the Event store : mapped and unmapped
  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
      mappedTrackCollection;

  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
      unmappedTrackCollection;

  // Map the whole event as one batch, which allows the mapper to process
  // the tracks concurrently
  std::vector<std::size_t> trackIds;
  std::vector<const Acts::RecordedMaterialTrack*> tracks;
  trackIds.reserve(mtrackCollection.size());
  tracks.reserve(mtrackCollection.size());
  for (const auto& [idTrack, mTrack] : mtrackCollection) {
    trackIds.push_back(idTrack);
    tracks.push_back(&mTrack);
  }

  auto results = m_cfg.materialMapper->mapMaterial(
      *m_mappingState, context.recoGeoContext, context.magFieldContext, tracks);

  for (std::size_t i = 0; i < results.size(); ++i) {
    auto& [mapped, unmapped] = results[i];
    mappedTrackCollection.try_emplace(mappedTrackCollection.end(), trackIds[i],
                                      std::move(mapped));
    unmappedTrackCollection.try_emplace(unmappedTrackCollection.end(),
                                        trackIds[i], std::move(unmapped));
  }

  // Write the mapped and unmapped material tracks to the output
  m_outputMappedMaterialTracks(context, std::move(mappedTrackCollection));
  m_outputUnmappedMaterialTracks(context, std::move(unmappedTrackCollection));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
