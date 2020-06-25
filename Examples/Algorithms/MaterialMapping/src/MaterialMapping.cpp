// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/MaterialMapping/MaterialMapping.hpp"

#include <iostream>
#include <stdexcept>

#include "ACTFW/Framework/WhiteBoard.hpp"

FW::MaterialMapping::MaterialMapping(const FW::MaterialMapping::Config& cnf,
                                     Acts::Logging::Level level)
    : FW::BareAlgorithm("MaterialMapping", level),
      m_cfg(cnf),
      m_mappingState(cnf.geoContext, cnf.magFieldContext),
      m_mappingStateVol(cnf.geoContext, cnf.magFieldContext) {
  if (!m_cfg.materialSurfaceMapper && !m_cfg.materialVolumeMapper) {
    throw std::invalid_argument("Missing material mapper");
  } else if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  ACTS_INFO("This algorithm requires inter-event information, "
            << "run in single-threaded mode!");

  if (m_cfg.materialSurfaceMapper) {
    // Generate and retrieve the central cache object
    m_mappingState = m_cfg.materialSurfaceMapper->createState(
        m_cfg.geoContext, m_cfg.magFieldContext, *m_cfg.trackingGeometry);
  }
  if (m_cfg.materialVolumeMapper) {
    // Generate and retrieve the central cache object
    m_mappingStateVol = m_cfg.materialVolumeMapper->createState(
        m_cfg.geoContext, m_cfg.magFieldContext, *m_cfg.trackingGeometry);
  }
}

FW::MaterialMapping::~MaterialMapping() {
  Acts::DetectorMaterialMaps detectorMaterial;

  if (m_cfg.materialSurfaceMapper && m_cfg.materialVolumeMapper) {
    // Finalize all the maps using the cached state
    m_cfg.materialSurfaceMapper->finalizeMaps(m_mappingState);
    m_cfg.materialVolumeMapper->finalizeMaps(m_mappingStateVol);
    // Loop over the state, and collect the maps for surfaces
    for (auto& [key, value] : m_mappingState.surfaceMaterial) {
      detectorMaterial.first.insert({key, std::move(value)});
    }
    // Loop over the state, and collect the maps for volumes
    for (auto& [key, value] : m_mappingStateVol.volumeMaterial) {
      detectorMaterial.second.insert({key, std::move(value)});
    }
  } else {
    if (m_cfg.materialSurfaceMapper) {
      // Finalize all the maps using the cached state
      m_cfg.materialSurfaceMapper->finalizeMaps(m_mappingState);
      // Loop over the state, and collect the maps for surfaces
      for (auto& [key, value] : m_mappingState.surfaceMaterial) {
        detectorMaterial.first.insert({key, std::move(value)});
      }
      // Loop over the state, and collect the maps for volumes
      for (auto& [key, value] : m_mappingState.volumeMaterial) {
        detectorMaterial.second.insert({key, std::move(value)});
      }
    }
    if (m_cfg.materialVolumeMapper) {
      // Finalize all the maps using the cached state
      m_cfg.materialVolumeMapper->finalizeMaps(m_mappingStateVol);
      // Loop over the state, and collect the maps for surfaces
      for (auto& [key, value] : m_mappingStateVol.surfaceMaterial) {
        detectorMaterial.first.insert({key, std::move(value)});
      }
      // Loop over the state, and collect the maps for volumes
      for (auto& [key, value] : m_mappingStateVol.volumeMaterial) {
        detectorMaterial.second.insert({key, std::move(value)});
      }
    }
  }
  // Loop over the available writers and write the maps
  for (auto& imw : m_cfg.materialWriters) {
    imw->writeMaterial(detectorMaterial);
  }
}

FW::ProcessCode FW::MaterialMapping::execute(
    const FW::AlgorithmContext& context) const {
  if (m_cfg.materialSurfaceMapper) {
    // Write to the collection to the EventStore
    std::vector<Acts::RecordedMaterialTrack> mtrackCollection =
        context.eventStore.get<std::vector<Acts::RecordedMaterialTrack>>(
            m_cfg.collection);

    // To make it work with the framework needs a lock guard
    auto mappingState =
        const_cast<Acts::SurfaceMaterialMapper::State*>(&m_mappingState);

    for (auto& mTrack : mtrackCollection) {
      // Map this one onto the geometry
      m_cfg.materialSurfaceMapper->mapMaterialTrack(*mappingState, mTrack);
    }

    context.eventStore.add(m_cfg.mappingMaterialCollection,
                           std::move(mtrackCollection));
  }
  if (m_cfg.materialVolumeMapper) {
    // Write to the collection to the EventStore
    std::vector<Acts::RecordedMaterialTrack> mtrackCollection =
        context.eventStore.get<std::vector<Acts::RecordedMaterialTrack>>(
            m_cfg.collection);

    // To make it work with the framework needs a lock guard
    auto mappingState =
        const_cast<Acts::VolumeMaterialMapper::State*>(&m_mappingStateVol);

    for (auto& mTrack : mtrackCollection) {
      // Map this one onto the geometry
      m_cfg.materialVolumeMapper->mapMaterialTrack(*mappingState, mTrack);
    }
  }
  return FW::ProcessCode::SUCCESS;
}
