// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/MaterialMapping/MaterialMapping.hpp"

#include "Acts/Material/AccumulatedMaterialSlab.hpp"
#include "Acts/Material/AccumulatedSurfaceMaterial.hpp"
#include "ActsExamples/MaterialMapping/IMaterialWriter.hpp"

#include <iostream>
#include <stdexcept>
#include <unordered_map>

namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples

ActsExamples::MaterialMapping::MaterialMapping(
    const ActsExamples::MaterialMapping::Config& cfg,
    Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("MaterialMapping", level), m_cfg(cfg) {
  if (!m_cfg.materialSurfaceMapper && !m_cfg.materialVolumeMapper) {
    throw std::invalid_argument("Missing material mapper");
  } else if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  m_inputMaterialTracks.initialize(m_cfg.collection);
  m_outputMaterialTracks.initialize(m_cfg.mappingMaterialCollection);

  ACTS_INFO("This algorithm requires inter-event information, "
            << "run in single-threaded mode!");

  if (m_cfg.materialSurfaceMapper) {
    // Generate and retrieve the central cache object
    m_mappingState = m_cfg.materialSurfaceMapper->createState(
        m_cfg.geoContext, m_cfg.magFieldContext);
  }
  if (m_cfg.materialVolumeMapper) {
    // Generate and retrieve the central cache object
    m_mappingStateVol = m_cfg.materialVolumeMapper->createState(
        m_cfg.geoContext, m_cfg.magFieldContext);
  }
}

ActsExamples::MaterialMapping::~MaterialMapping() {
  Acts::DetectorMaterialMaps detectorMaterial;

  if (m_cfg.materialSurfaceMapper != nullptr &&
      m_cfg.materialVolumeMapper != nullptr) {
    // Finalize all the maps using the cached state
    auto surfaceDetectorMaterial =
        m_cfg.materialSurfaceMapper->finalizeMaps(*m_mappingState);
    auto volumeDetectorMaterial =
        m_cfg.materialVolumeMapper->finalizeMaps(*m_mappingStateVol);
    // Loop over the state, and collect the maps for surfaces
    for (auto& [key, value] : surfaceDetectorMaterial.first) {
      detectorMaterial.first.insert({key, std::move(value)});
    }
    // Loop over the state, and collect the maps for volumes
    for (auto& [key, value] : volumeDetectorMaterial.second) {
      detectorMaterial.second.insert({key, std::move(value)});
    }
  } else if (m_cfg.materialSurfaceMapper != nullptr) {
    detectorMaterial =
        m_cfg.materialSurfaceMapper->finalizeMaps(*m_mappingState);
  } else if (m_cfg.materialVolumeMapper != nullptr) {
    detectorMaterial =
        m_cfg.materialVolumeMapper->finalizeMaps(*m_mappingStateVol);
  }
  // Loop over the available writers and write the maps
  for (auto& imw : m_cfg.materialWriters) {
    imw->writeMaterial(detectorMaterial);
  }
}

ActsExamples::ProcessCode ActsExamples::MaterialMapping::execute(
    const ActsExamples::AlgorithmContext& context) const {
  // Take the collection from the EventStore
  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
      mtrackCollection = m_inputMaterialTracks(context);

  if (m_cfg.materialSurfaceMapper) {
    // To make it work with the framework needs a lock guard
    for (auto& [idTrack, mTrack] : mtrackCollection) {
      // Map this one onto the geometry
      m_cfg.materialSurfaceMapper->mapMaterialTrack(*m_mappingState, mTrack);
    }
  }
  if (m_cfg.materialVolumeMapper) {
    // To make it work with the framework needs a lock guard
    for (auto& [idTrack, mTrack] : mtrackCollection) {
      // Map this one onto the geometry
      m_cfg.materialVolumeMapper->mapMaterialTrack(*m_mappingStateVol, mTrack);
    }
  }
  // Write take the collection to the EventStore
  m_outputMaterialTracks(context, std::move(mtrackCollection));
  return ActsExamples::ProcessCode::SUCCESS;
}

std::vector<std::pair<double, int>>
ActsExamples::MaterialMapping::scoringParameters(uint64_t surfaceID) {
  std::vector<std::pair<double, int>> scoringParameters;

  if (m_cfg.materialSurfaceMapper) {
    Acts::SurfaceMaterialMapper::State* smState =
        static_cast<Acts::SurfaceMaterialMapper::State*>(m_mappingState.get());

    auto surfaceAccumulatedMaterial =
        smState->accumulatedMaterial.find(Acts::GeometryIdentifier(surfaceID));

    if (surfaceAccumulatedMaterial != smState->accumulatedMaterial.end()) {
      auto matrixMaterial =
          surfaceAccumulatedMaterial->second.accumulatedMaterial();
      for (const auto& vectorMaterial : matrixMaterial) {
        for (const auto& AccumulatedMaterial : vectorMaterial) {
          auto totalVariance = AccumulatedMaterial.totalVariance();
          scoringParameters.push_back(
              {totalVariance.first, totalVariance.second});
        }
      }
    }
  }
  return scoringParameters;
}
