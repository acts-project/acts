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
                           << "run in single-threaded mode!");

  if (m_cfg.materialMapper == nullptr) {
    throw std::invalid_argument("Missing material mapper");
  }
  // Create the state object
  m_mappingState = m_cfg.materialMapper->createState(m_cfg.geoContext);
}

MaterialMapping::~MaterialMapping() {
  Acts::TrackingGeometryMaterial detectorMaterial =
      m_cfg.materialMapper->finalizeMaps(*m_mappingState, m_cfg.geoContext);
  // Loop over the available writers and write the maps
  for (auto& imw : m_cfg.materialWriters) {
    imw->writeMaterial(detectorMaterial);
  }
}

ProcessCode MaterialMapping::execute(const AlgorithmContext& context) const {
  // Take the collection from the EventStore: input collection
  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
      mtrackCollection = m_inputMaterialTracks(context);

  // Write the output collections to the Event store : mapped and unmapped
  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
      mappedTrackCollection;

  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
      unmappedTrackCollection;

  for (const auto& [idTrack, mTrack] : mtrackCollection) {
    auto [mapped, unmapped] = m_cfg.materialMapper->mapMaterial(
        *m_mappingState, context.geoContext, context.magFieldContext, mTrack);

    mappedTrackCollection.try_emplace(mappedTrackCollection.end(), idTrack,
                                      mapped);
    unmappedTrackCollection.try_emplace(unmappedTrackCollection.end(), idTrack,
                                        unmapped);
  }

  // Write the mapped and unmapped material tracks to the output
  m_outputMappedMaterialTracks(context, std::move(mappedTrackCollection));
  m_outputUnmappedMaterialTracks(context, std::move(unmappedTrackCollection));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
