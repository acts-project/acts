// This file is part of the Acts project.
//
// Copyright (C) 2017-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/MaterialMapping/MaterialMapping.hpp"

#include "Acts/Material/AccumulatedMaterialSlab.hpp"
#include "Acts/Material/AccumulatedSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "ActsExamples/MaterialMapping/IMaterialWriter.hpp"

#include <stdexcept>
#include <unordered_map>

namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples

ActsExamples::MaterialMapping::MaterialMapping(
    const ActsExamples::MaterialMapping::Config& cfg,
    Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("MaterialMapping", level), m_cfg(cfg) {
  // if (m_cfg.materialMappers.empty()) {
    // throw std::invalid_argument(
        // "At least one material mapper must to be defined.");
  // } else if (m_cfg.trackingGeometry == nullptr and m_cfg.detector == nullptr) {
    // throw std::invalid_argument(
        // "Either TrackingGeometry or Detector must to be defined.");
  // }

  if (m_cfg.trackingGeometry == nullptr and m_cfg.detector == nullptr) {
    throw std::invalid_argument(
        "Either TrackingGeometry or Detector must to be defined.");
  }

  if (m_cfg.materialSurfaceMapper != nullptr) {
    m_cfg.materialMappers.push_back(m_cfg.materialSurfaceMapper);
  }
  if (m_cfg.materialVolumeMapper != nullptr) {
    m_cfg.materialMappers.push_back(m_cfg.materialVolumeMapper);
  }

  if (m_cfg.trackingGeometry != nullptr) {
    ACTS_DEBUG("Material mapping for TrackingGeometry is enabled.");
    for (const auto& mapper : m_cfg.materialMappers) {
      m_mappingStates.push_back(mapper->createState(
          m_cfg.geoContext, m_cfg.magFieldContext, *m_cfg.trackingGeometry));
    }

  } else if (m_cfg.detector != nullptr) {
    ACTS_DEBUG("Material mapping for Detector is enabled.");
    for (const auto& mapper : m_cfg.materialMappers) {
      m_mappingStates.push_back(mapper->createState(
          m_cfg.geoContext, m_cfg.magFieldContext, *m_cfg.detector));
    }
  }

  m_inputMaterialTracks.initialize(m_cfg.collection);
  m_mappedMaterialTracks.initialize(m_cfg.mappedCollection);
  m_unmappedMaterialTracks.initialize(m_cfg.unmappedCollection);
}

ActsExamples::MaterialMapping::~MaterialMapping() {
  // Prepare to write out the maps
  Acts::DetectorMaterialMaps detectorMaterial;
  // Finalize all the mappers and collect the maps
  for (const auto [im, mapper] : Acts::enumerate(m_cfg.materialMappers)) {
    Acts::MaterialMappingState& mState = *m_mappingStates[im].get();
    auto mappingResult = mapper->finalizeMaps(mState);
    // The surface maps - check for duplicates
    for (auto& [id, surfaceMaterial] : mappingResult.surfaceMaterial) {
      if (detectorMaterial.first.find(id) != detectorMaterial.first.end()) {
        ACTS_ERROR("Surface material map already exists for this identifier.");
      }
      detectorMaterial.first.insert({id, std::move(surfaceMaterial)});
    }
    // The volume maps - check for duplicates
    for (auto& [id, volumeMaterial] : mappingResult.volumeMaterial) {
      if (detectorMaterial.second.find(id) != detectorMaterial.second.end()) {
        ACTS_ERROR("Volume material map already exists for this identifier.");
      }
      detectorMaterial.second.insert({id, std::move(volumeMaterial)});
    }
  }
  // Loop over the available writers and write the maps
  for (auto& writer : m_cfg.materialWriters) {
    writer->writeMaterial(detectorMaterial);
  }
}

ActsExamples::ProcessCode ActsExamples::MaterialMapping::execute(
    const ActsExamples::AlgorithmContext& context) const {
  // Take the input collection from the EventStore
  std::unordered_map<size_t, Acts::RecordedMaterialTrack> inputMaterialTracks =
      m_inputMaterialTracks(context);

  // Write the output collection to the EventStore
  std::unordered_map<size_t, Acts::RecordedMaterialTrack> mappedMaterialTracks;
  std::unordered_map<size_t, Acts::RecordedMaterialTrack>
      unmappedMaterialTracks;

  for (auto [idTrack, mTrack] : inputMaterialTracks) {
    // The processed part of the track
    Acts::RecordedMaterialTrack mappedTrack{mTrack.first,
                                            Acts::RecordedMaterial{}};
    Acts::RecordedMaterialTrack unmappedTrack{mTrack.first,
                                              Acts::RecordedMaterial{}};
    // Map this one onto the geometry
    for (const auto [im, mapper] : Acts::enumerate(m_cfg.materialMappers)) {
      Acts::MaterialMappingState& mState =
          const_cast<Acts::MaterialMappingState&>(*(m_mappingStates[im].get()));
      auto [mapped, unmapped] = mapper->mapMaterialTrack(mState, mTrack);
      // Run further with the unmapped part
      mTrack = unmapped;
      // Keep track of the mapped part
      mappedTrack.second.materialInX0 += mapped.second.materialInX0;
      mappedTrack.second.materialInL0 += mapped.second.materialInL0;
      mappedTrack.second.materialInteractions.insert(
          mappedTrack.second.materialInteractions.end(),
          mapped.second.materialInteractions.begin(),
          mapped.second.materialInteractions.end());
      // The unmapped is handed through to the end
      unmappedTrack = unmapped;
    }
    mappedMaterialTracks.insert({idTrack, mappedTrack});
    unmappedMaterialTracks.insert({idTrack, unmappedTrack});
  }
  // Write the result to the event store
  m_mappedMaterialTracks(context, std::move(mappedMaterialTracks));
  m_unmappedMaterialTracks(context, std::move(unmappedMaterialTracks));
  // Return success
  return ActsExamples::ProcessCode::SUCCESS;
}